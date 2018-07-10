
## Instructions for Graphing JaBba Outputs using gGnome/gGnome.js

1. Make sure on the cluster you have both gGnome and JaBba installed. The instructions for setting up these packages can be found at the following websites:

    [JaBba Github Page](https://github.com/mskilab/JaBba)

    [gGnome Github Page](https://github.com/mskilab/gGnome)

2. Make sure gGnome.js is installed on your local machine. It is recommended that you clone this project into your computer's lowest directory, i.e. opening terminal and cloning without changing directory. The instructions for setting up this package can be found here:

   [gGnome.js Github Page](https://github.com/mskilab/gGnome.js)

3. Once these three packages are installed, run JaBba on your junctions and coverage file. This can be run from an `R` session or using the JaBba executable. Be sure to specify the output directory as this is where your input `.rds` file will be for the following steps.

  ```R
  ## If running from an R session use this code, if running executable don't worry about this  
  JaBbA(junctions = "Path to junctions file",  
        coverage = "Path to coverage file",  
        outdir = "Directory to dump into")  
  ```

4. Load gGnome into an open `R` session. Build a `gGnome` object from your `JaBbA` output file `jabba.simple.gg.rds`. You may also build a `gGnome` object from your output file `jabba.simple.rd` using `gread` as specified below.

  ```R
  library(gGnome)  

  ## Use if loading jabba.simple.gg.rds  
  gg = gGraph$new(jabba = "Path to jabba.simple.gg.rds")  
  
  ## Use if loading jabba.simple.rds  
  gg = gread("Path to jabba.simple.rds")
  ```

5. Now, convert your `gGraph` object to a `.json` file by running the following code:

  ```R
  ## To get a standard graph  
  gg$gg2js(filename = "yourfilename.js")  

  ## If you would like to graph your gGraph without copy number  
  gg$gg2js(filename = "yourfilename.js", no.y = TRUE)  
  ```

6. Step 5 should have generated a file named `"yourfilename.js"` in your current working directory. Navigate to this file in the cluster using Finder on your local machine (make sure you mounted `mount.nygc.sh`). Move this file into the directory `gGnome.js/json`. This directory should be located within your home on your local machine.

7. Start the gGnome.js project by opening a new terminal window and running the start file. This should open a new window in your Chrome browser called `http://localhost:8080/index.html`. If this does not open within a minute, click on the hyperlink in this step.

  ```bash
  $ cd gGnome.js/  
  $ ./start.sh  
  ```

8. From within the open window for gGnome.js, select your file from the drag down list of files located at the top of your screen. You now have loaded a JaBbA output into the gGnome.js visualization! Go cure cancer!
