## Singularity RStudio and Shiny Server

Image based on the Ubuntu 20.04 (R 4.2.2), see `./singularity/rstudio_server/recipe-4-2-2.def` for details

### First step (get the image)

Download the image at <https://cloud.sylabs.io/library/michael.rade/r-world/rstudio-server> or pull it:

``` sh
$ singularity pull --arch amd64 library://michael.rade/r-world/rstudio-server:1.0.0
```

### Second step (config and run)

The following adjustments must be made in the `run.sh` file:

-   Set a port (variable `RSTUDIO_PORT`)
-   Set a path for the log files (variable `TMPDIR`)
-   The variable `RSTUDIO_PASSWORD` can be set to a password of your choice. Default is "password"
-   For parameter `auth-pam-helper` and `rsession-config-file` in run.sh: Set path to this folder.

Furthermore, set your R library path in `./singularity/rsession.conf`. I recommend to create a new folder. It can lead to conflicts when folders with libraries of another R version are used.

Now run the following :

``` sh
$ bash run.sh run rstudio-server:1.0.0.sif
```

You can reach RStudio Server via your webbrowser (e.g. `localhost:8072/auth-sign-in`). Use the port that you have specified in the variable `RSTUDIO_PORT`. You can log in with your current user name and password you set in `RSTUDIO_PASSWORD`.

### Notes

You can also work on the command line

``` sh
# If you need to bind something that is not your user's home directory
$ export SINGULARITY_BINDPATH="/mnt/example/"
# Define path to R libraries
$ export R_LIBS_USER="/home//michael.rade/R/x86_64-pc-linux-gnu-library/4.0.2.RStudioServer/"

$ singularity shell rstudio-server:1.0.0.sif
```

You can build a Singularity image with:

``` sh
$ sudo singularity build rstudio-server:1.0.0.sif recipe-4-2-2.def
```

You can put more parameter to run.sh. To see all possible parameter use:

``` sh
$ singularity run --app rserver rstudio-server_1.0.0.sif --help
```
