# Notes and thoughts during package development

A place to track my notes and thoughts as I work on converting the scripts to an R package. 

## Current gameplan

In thinking about the the best way to approach this, I think I've settled on a relatively minimal approach to start. I'm going to do relatively little code re-writing to start. First, I want to just port the long script of functions we currently have into a new package structure, and to document, as much as I can, the "user-facing" functions in the pipeline. Along the way, I'd like to improve the error messages for these a little bit, tidy up some small things (fewer print statements, make dependencies clearer).

To do this, I'll try to separate out the functionality of each of the "main" user-facing functions (that is, the ones that are called in the example pipelines) into a single .R file. That file will define the main function, and all the "subfunctions"" it uses. For now, I think I'll export all functions, but maybe in the future we won't need to do that. I'll preface the name of each .R file with "old", to mark it as the "raw" input. 

Then, maybe people can use the package for a little while. Find some bugs or issues, talk about what they like and don't like, etc. Then, can attempt to do some re-writing, if needed, based on that user feedback.

One issue that hangs over all this is testing. I think testing would be great to do, but I may hold off on it a little bit for now. I think the code would maybe need to be re-written to make it more testable. And maybe I just need to get it all in my head a little more before I think about the best way to test. 


### Tidbits 

Working from the Kintler DIA data, have already run into an issue with the data import: missing one one of the necessary input columns. Error message doesn't say which, should change that. In checking, seems to be the "exclusivity" column. Stephanie suggests we maybe don't use that column anyway? For now, just added that column (blank) to the Kintler DIA data for testing. Upon doing that, the extract_data function works. 


#### 2022-04-04

Working through documentation and very minor code editing (e.g., removing extra print object declarations for printing behind semicolons) of extract_data subfunctions for DIA data. 
