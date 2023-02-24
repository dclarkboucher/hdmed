## Test Enviroments
* Local Windows install, R 4.2.2
* Windows devel, release (4.2.2), and oldrelease (4.1.3)
* OS X release

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:

  None of these words are misspelled - they are either Latin abbreviations
  (e.g., "et al.") or acronyms, all of which are defined in the DESCRIPTION file text
  
## Downstream dependencies
There are no downstream dependencies of this package

## Resubmission
This is a resubmission. In this version, I have:

* Replaced a call to citEntry() in the CITATION file to bibentry()

* Removed a call to personList() in the CITATION
