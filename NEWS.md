# drugdevelopR 1.0.2

Improvements that will work with the Shiny apps
* Adding trace to all `optimal_*()` functions.

# drugdevelopR 1.0.1

Incorporating suggestions during CRAN submission: 
* Replacing `\dontrun{}` by `\donttest{}` in all examples that can be run in theory, but are too long to be run in practice.
* Replacing non-suppressable output to the console by a custom print function for output objects and customizable progress bars using progressr.
* We removed function code setting the seed to a specific number in R/functions_multiple_normal.R; R/functions_multiple_tte.R.

# drugdevelopR 1.0.0

* Complete validated version with full documentation of all user-relevant functions.

# drugdevelopR 0.1.0

* Added a `NEWS.md` file to track changes to the package.
