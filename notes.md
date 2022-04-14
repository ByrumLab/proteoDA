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

#### 2022-04-05

One architecture thing I'm noticing that doesn't seem ideal: the annotation and individual intensity data get split up into separate data frames, with the protein ID data in the intensity data getting put into the rownames. Looks like these get checked again in the future to make sure they still match up and order hasn't changed for some reason: the first steps of the run_limma_analysis function do this. But, maybe better to just not break these apart in the first place? I guess some of that will depend on what happens down the line: not sure if some operations later on need to operate on just the matrix of sample data without the protein annotation, making it somewhat easier if these are split. In any case, looks like down the line it is just a check: there isn't currently any code to fix this is it isn't true. 

Related to that, a lot of information seems to get put into rownames, instead of into a column. Need to double-check if there's some performance or particular reason for that, but my personal preference would be to put that into back in columns, so that it can be accessed in the normal way ($, [], etc) instead of having to use rownames(). My intuition is that we shouldn't really use rownames at all: just use a normal named character column. That also provides more detail: the column name will provide more info (hopefully) about what the value actually is: e.g., exactly what type of name or protein ID it is. 

More architecture ideas: I think the current overall idea, with high-level wrapper functions that call a bunch of subfunctions, is the right idea in general. Right now, though, this is a little tough to reason about, as the high-level functions both (1) call subfunctions, and (2) do other data processing steps just within them. I think reasoning about the control flow of information through the function would be easier if the toplevel wrapper was really just that, and if the processing steps that are currently done in the function body were moved into subfunctions. This would also (hopefully) make customization and debugging easier: if you ran into a bug running the wrapper function, you could run each of the subfunctions in turn, but wouldn't need to worry about any extra code. 

Another general code issue I'm finding: too-short, nondescriptive names. The output list structures also seem redundant. Often, some of the slots in the output list are just summaries of info in other parts of the list: e.g., the number of proteins is in one slot, but you could also get that just by looking at the number of rows in the protein data frame. Not sure how I feel about this redundancy: could sometimes be helpful when the list structures are unclear: these elements are better named than some of the other ones. 

Another general issue to think about: how and whether to separate pipelines. As it stands now, most functions have a few if statement blocks to switch between the different tasks for each pipeline. Not sure the best way to architect this overall. For the user-facing function, nice to just have one function and put in pipeline as an argument. All the if statements are a little tough, though. 

I do (possibly) like the idea of using lists with some slots for pipe and enrichment, such that it becomes a sort of weak type to deal with multiple dispatch and making it easier to pass the pipe and enrichment options down through all the subfunctions. 


#### 2022-04-14

A re-architecture idea: disentangle and clarify the tasks of the first few functions in the pipeline. Right now, they all step on each other in unclear ways. Extract_data does data extraction, but also does some processing (e.g., remove contaminants, which removes rows from the data file). subset_targets filters out samples from the targets file. Then, process_data does both data filtering (removes rows from data) and sample filtering (removes columns from data), while taking on normalization as well. 

I think I would rework things. An import data function, that pulls in the Maxquant data but does little else. A simple function that imports the metadata as a dataframe. Then, a function that processes the metadata into a targets file (where you have the option to remove samples if you want). Then, a function that takes in the Maxquant data and the targets dataframe or file and does all the cross-checking against each other, sample filtering, validation, etc. Then, a row/protein filtering step. Then, a normalization step. Still need to think on this more (not as sure about the combination step stuff). But I do think some conceptual separation of all these early functions would be helpful. 

