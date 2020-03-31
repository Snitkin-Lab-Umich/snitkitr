# snitkitr
## This R package contains functions that are useful across projects/users. 

## Getting the snitkitr package onto Great Lakes (The University of Michigan HPC cluster)

`$ R # Start R`  
`> devtools::install_github("Snitkin-Lab-Umich/snitkitr", force = TRUE)`  
`> library(snitkitr)`  
  
You can now quit R. In your Rscripts call `library(snitkitr)` at the beginning or use `snitkitr::INSERT_FUNCTION_NAME().`. 

## Keeping the snitkitr package up-to-date  
Setting `force = TRUE` in the above command requires R to get the most up-to-date version of snitkitr from github.   

## Adding or editing a snitkitr package 
To add a function or edit a function in the package clone the github repo to your local computer. To clone a repo you need to have a github account. Once you have an account click on the "Clone or download" button on this page (https://github.com/Snitkin-Lab-Umich/snitkitr). Copy the text "https://github.com/Snitkin-Lab-Umich/snitkitr.git" that appears when you click the button. 

On your local computer open a terminal and head to the location where you want your local clone of the repository to live. Once you are in the desired location:   
`$ git clone https://github.com/Snitkin-Lab-Umich/snitkitr.git`  
`$ cd snitkitr/R`  
  
Now you're ready to add a function. If your function is related to currently available functions add it to the end of the appropriate file. If your function is novel start a new file.   

Every time you go to the snitkitr/ directory run:  
`$ git pull`   
This command will make sure you have any changes other users have made to the package on your local copy.  
  
Every time you're done writing code run:   
`$ git status` # To see which files you've edited   
`$ git add <file name>`  # to add multiple files use `$ git add .`  
`$ git commit -m "In present tense, add a ~one sentence description of the changes you made"` # This is a record to explain to yourself and others what code you wrote   
`$ git push` # This shares your changes to github so that other users can grab them  
  
## Writing code that's easy for others to use
When writing functions for the whole lab to use please make your code easy to read. There are several steps you can take to make the code accessible to other users:  
1. Write informative function and variable names. 
2. Use the R style guide (http://adv-r.had.co.nz/Style.html). 
3. Comment your code.
4. Check that inputs to functions are of the type and size expected.
5. Writing warnings or error messages if inputs or results are incorrect. 
6. Write unit tests. 
7. Provide an example. 
8. Break long functions up into more descriptive subfunctions. 

## Specifically, please comment your function using the R package function documentation style called Roxygen2. 
1. In RStudio put your cursor on the name of the function you've written. 
2. Code > Insert Roxygen Skeleton
3. Update the title to be a ~one line phrase describing the function. 
4. On the next line add @description and write ~one paragraph describing the function in more detail. 
5. For each input to the function describe it in the @param section. 
6. Describe whatever gets returned by the function in @return. 
7. No need to write anything after @export, but make sure it stays in the documentation. This command allows users to access this function. 
8. Write examples using the function, if applicable. 

## Adding a unit test for a function 
Unit tests are standarized ways to check that a function behaves as expected. 
To write a unit test for a function you've written open snitkitr as a .Rproj. 
Edit or create a test-R_FILE_NAME.R file in snitkitr/tests/testhat/  

For an example see test-greatlakes_resource_usage.R:  
```test_that("ConvertTimeToHours returns expected times", {
  time <- "10:00:00"
  out <- ConvertTimeToHours(time)
  expect_equal(10, out)
  
  time <- "10:30:00"
  out <- ConvertTimeToHours(time)
  expect_equal(10.5, out)
  
  time <- "1-10:00:00"
  out <- ConvertTimeToHours(time)
  expect_equal(34, out)
})

test_that("ConvertTimeToHours gives error with incorrect input", {
  time <- 10
  expect_error(ConvertTimeToHours(time))
  
  time <- "10:30"
  expect_error(ConvertTimeToHours(time))
})
```

The `expect_equal()` function tests that the two inputs are equal; if they are not the unit test would fail.
The `expect_error()` function tests that what is called within gives an error or stop. 

There are other `expect_()` functions like `expect_true()`, `expect_false()`, etc... 

## Running the unit tests
To check if all of the unit tests are working, which means all of the package functions are working as expected click: 
Build > More > Test Package
Wait to see if all of the tests pass. If some fail you'll need to update your function code. 


