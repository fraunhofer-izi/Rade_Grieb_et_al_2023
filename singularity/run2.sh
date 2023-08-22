#!/bin/sh

module load Singularity # or install singularity: https://docs.sylabs.io/guides/latest/user-guide/quick_start.html
export SINGULARITY_BINDPATH="/mnt/ribolution/,/mnt/workdata1/" # You can remove it when there is nothing to bind

export RSTUDIO_PORT=8062

# we don't want to clutter the cluster directories, so each user sets up his own
# rstudio-server-tmp for now
export TMPDIR="/homes/olymp/$USER/RStudio-Server-Files/"

echo "Path where RStudio Server will write run-time state, logs, etc.:"
echo $TMPDIR

mkdir -p "$TMPDIR/tmp/rstudio-server"
uuidgen > "$TMPDIR/tmp/rstudio-server/secure-cookie-key"
# the secure cookie should of course only be accessible by the owner
chmod 0600 "$TMPDIR/tmp/rstudio-server/secure-cookie-key"

mkdir -p "$TMPDIR/var/lib"
mkdir -p "$TMPDIR/var/run"

RSTUDIO_PASSWORD="password" singularity "$1" \
  --bind="$TMPDIR/var/lib:/var/lib/rstudio-server" \
  --bind="$TMPDIR/var/run:/var/run/rstudio-server" \
  --bind="$TMPDIR/tmp:/tmp" \
  "$2" \
  --server-user $USER \
  --www-port $RSTUDIO_PORT \
  --auth-none 0 \
  --auth-pam-helper ~/projects/2020-imSAVAR/singularity/rstudio_server/rstudio_auth.sh \
  --rsession-config-file ~/projects/2020-imSAVAR/singularity/rstudio_server/rsession.conf \
  --auth-timeout-minutes=0 \
  --auth-stay-signed-in-days=7
