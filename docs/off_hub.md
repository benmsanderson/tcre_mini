# Running the WIEMIP analysis off the JupyterHub

The JupyterHub is not needed for any step. It provides two things, both of which another machine can
reproduce:

- **The data.** `/srv/wiemip` on the hub is the Wasabi bucket `wiemip` mounted read-only with s3fs
  (endpoint `https://s3.us-east-2.wasabisys.com`, region `us-east-2`, path-style requests).
- **The registry.** `wiemip_registry` is installed from `WIEMIP/wiemip-data-registry` (private;
  `scripts/check_wiemip_access.py` prints the commit in use). It finds the data through two environment
  variables: `WIEMIP_DATA_ROOT` (default `/mnt/wiemip`) and `WIEMIP_CSV_PATH` (default `/srv/wiemip-csv`,
  a cache of global sums that the registry fills on first read).

Only `scripts/wiemip_extract.py` reads WIEMIP files. Everything after it (`wiemip_values.py`,
`wiemip_power.py`, `wiemip_paired_spread.py`, `wiemip_paired_tests.py`) reads the CSV files in
`data/wiemip/` and `results/wiemip/`, so the long PyMC runs need only this repository.

WIEMIP output is unpublished. Keep the local copy and the bucket keys private, and share derived data
only with WIEMIP authors.

## 1. Environment

```bash
git clone https://github.com/benmsanderson/tcre_mini && cd tcre_mini   # private; needs access
git submodule update --init external/TCRE_Assessment
conda env create -f environment.yml && conda activate tcre-evidence
python -m pytest -q
```

For the extraction step only, install the registry into the same environment (it needs access to the
WIEMIP GitHub organisation; use the commit `check_wiemip_access.py` reports on the hub to match it):

```bash
pip install "git+https://github.com/WIEMIP/wiemip-data-registry.git@<commit>"
```

## 2. Bucket credentials

Set these up yourself; do not paste keys into a chat or commit them.

- For s3fs: a file `~/.passwd-s3fs` containing `ACCESS_KEY:SECRET_KEY`, with `chmod 600 ~/.passwd-s3fs`.
- For the AWS CLI or rclone: `aws configure --profile wiemip` (any region; the endpoint is given below).

## 3. Make the data readable

Either mount the bucket (Linux, s3fs; on macOS `rclone mount` is simpler):

```bash
mkdir -p ~/wiemip
s3fs wiemip ~/wiemip -o ro,use_path_request_style,passwd_file=$HOME/.passwd-s3fs \
     -o url=https://s3.us-east-2.wasabisys.com,endpoint=us-east-2
```

or copy the 1pctCO2 part once, which avoids network failures during a run (check the size first):

```bash
E=https://s3.us-east-2.wasabisys.com
aws s3 ls s3://wiemip/1pctCO2/ --recursive --summarize --human-readable --endpoint-url $E --profile wiemip | tail -2
aws s3 sync s3://wiemip/1pctCO2/ ~/wiemip/1pctCO2/ --endpoint-url $E --profile wiemip
```

Then point the registry at it:

```bash
export WIEMIP_DATA_ROOT=~/wiemip WIEMIP_CSV_PATH=~/wiemip-csv
```

A cold cache makes the first extraction slow, because every global sum is computed from the gridded
files. Copying the hub's `/srv/wiemip-csv` (about 170 MB) to `$WIEMIP_CSV_PATH` avoids that, provided the
registry commit is the same.

## 4. Check, extract, fit

```bash
python scripts/check_wiemip_access.py --read     # paths, inputs, one probe per model; exit 1 on failure
python scripts/wiemip_extract.py                 # writes data/wiemip/
python scripts/wiemip_values.py && python scripts/wiemip_power.py
python scripts/wiemip_paired_spread.py --ensemble wiemip --jobs <n>
python scripts/wiemip_paired_tests.py synthetic --jobs <n>
python scripts/wiemip_paired_tests.py loo --jobs <n>
```

`wiemip_paired_spread.py` and `wiemip_paired_tests.py` take minutes to about 1.5 hours on 16 cores. If the
machine can drop your session, start them detached and log to a file, for example
`setsid nohup python scripts/wiemip_paired_tests.py synthetic --jobs 16 > run.log 2>&1 &`.
`wiemip_paired_tests.py` writes its results after each quantity; `wiemip_paired_spread.py` writes at the end.

## 5. Without the bucket

If only the fits are needed, skip steps 2 and 3 and the first two commands of step 4: the committed
`data/wiemip/` and `results/wiemip/` files are enough to run `wiemip_values.py` onwards.
