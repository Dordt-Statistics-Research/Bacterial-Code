README - documentation for Bacterial code as of 9/12/16 (author Craig Disselkoen)

-------------------------------------------------------
---------------------*** CODE ***----------------------
-------------------------------------------------------


-------------------------------------------------------
The following R files constitute the "core" of the Bacterial project and the public code release for Bayesian 2.0.

  aijContext.R : Functions for constructing and dealing with 'aijContext' objects (using R's S3 object system).
    
    R doesn't really have a way to mark 'public' and 'private' fields/methods, but the intention here is 
    that all fields / members are private, and (outside of the aijContext.R file) aijContext objects should
    only be accessed using methods.  Methods intended to be public and methods intended to be private
    are clearly marked as such in the file, and additionally, the public methods have signatures and 
    documentation provided at the top of this file.  
  
  gibbs_sampler_multivariate.R : Code implementing the Gibbs sampler itself, called in aijContext.R. 
    
    This code was originally written by Brian Greco, and I don't actually understand what it's doing all
    the time.  I have optimized it a little, though. 
  
  Aij_generation.R : Implementations of our aij-generating methods, as well as competing methods we benchmark against.
    
    Really, this file could be viewed as an extension of aijContext.R, as it defines even more methods
    on aijContext objects.  However, it seemed good to isolate these key methods in their own file.
    Also, it keeps the aijContext.R file from becoming insanely long.

-------------------------------------------------------
The following R files implement the handling of our own input data files.
They also implement automatic caching functionality, described more completely below.

  wrAijContext.R : Functions for constructing and dealing with 'wrAijContext' ("wrapped aijContext") objects.
    
    wrAijContext is intended to be a subclass of aijContext.  I didn't actually use R's S3 inheritance system
    because (a) I just didn't, and (b) I found it valuable to have to explicitly list all the methods which were
    being inherited, so that it is clear which methods are explicitly being inherited without modifications
    and which are being modified/wrapped.  This way, you can't forget a method and 'inherit' it without realizing.
    
    Like aijContext objects, wrAijContext objects are intended to be treated as having private fields/members,
    and only be accessed using the methods listed in wrAijContext.R.  Note that wrAijContext.R adjusts the 
    implementation of many of the fundamental accessor methods from aijContext.  This makes it doubly important
    that aijContext methods (including even the public methods in aijContext.R) use only the accessor methods
    to access aijContext objects.  As a result, all of the public aijContext methods accept aijContext or 
    wrAijContext objects equally, and will work correctly in either case.  
    
    wrAijContext objects are extremely useful upgrades to aijContext objects.  First of all, there is a 
    constructor for wrAijContext (get.wrapped.aij.context.from.ds) which simply takes a string describing
    a data source, and creates a wrAijContext automatically.  (There are other constructors as well.)
    More importantly, wrAijContext objects have automatic caching.  This results in huge performance gains
    vs. aijContext objects.  Basically, whenever you make a computationally-intensive call to a method on a
    wrAijContext object, the result is saved away (in the Rsave/ subdirectory) before it is returned.
    Then, in the future, if you make the same call with the same wrAijContext object, the result is simply
    loaded from file and returned to you, instead of recomputed.  More info, explanation, and documentation
    is available at the top of the wrAijContext.R file. 

  Aij_generation_input_data.R : Functions for dealing with specifically our input data files
  
    We have a lot of input data, and a lot of different varieties.  This contains the code to (basically)
    parse each possible file and return the contents in a uniform format.  These functions are used by 
    wrAijContext objects, and are also useful on their own. 
    
  genenameMaps.R : Supplementary script which generates maps for converting between various gene name formats.

-------------------------------------------------------
The following R files are mostly small, useful utility functions for working with the project.

  expnameMaps.R : Like genenameMaps.R above, but for experiment names.
  
    Because experiment-name mappings are many-to-one rather than one-to-one (like gene name mappings),
    dealing with experiment names is a little more complicated; hence, expnameMaps.R provides methods
    for looking up names in the mapping, rather than just giving you a data frame like genenameMaps.R. 
  
  overallmap.Rdata : Not a script, but a saved data frame.
  
    Basically contains a better genename mapping that I put together at some point.  It's the same as
    would be obtained from genenameMaps.R except for the common names.  We have two different sources of
    mappings to common names, which sometimes disagree; this overallmap.Rdata includes both mappings and 
    a column simply called "common" which includes the consensus common name for each gene, or in case of 
    disagreement, both common names.  Some of the scripts here use overallmap.Rdata. 
  
  Aij_utility_funcs.R : Useful methods for making your life easier, and creating basic graphs
  
    A handful of useful functions for common and repetitive tasks; and code for generating the 
    most-often-used graphs of various data.  Notably, these functions do not depend on valContext.R
    (below), and will accept either aijContext objects or wrAijContext objects. 

  craig_utility_funcs.R : Useful functions for general use in R coding
  
    A few of my personal favorite functions for coding in R.  None of these are bacterial-specific at all,
    so they can be used in any R coding project.  Included here because Aij_utility_funcs.R references a 
    few of them. 
  
-------------------------------------------------------
The following R files implement functions for examining the data produced by the above scripts,
and validating the performance of various methods vs. each other. 
  
  valContext.R : Functions for constructing and dealing with 'valContext' ("validation context") objects.
  
    valContext is intended to be a subclass of wrAijContext, in the same way that wrAijContext is a subclass
    of aijContext.  The key difference between valContext and wrAijContext is that a valContext requires
    "gold standard" answers regarding gene activity states (or some data meant to be used as a gold standard).
    If you understand wrAijContext.R, you shouldn't have any problems at all understanding this file. 
  
  Aij_validation.R : Lots of functions for outputting various statistics and graphs
  
    These functions mostly operate on valContext objects and produce various statistics and graphs related
    to validation.  I've rarely ever deleted anything from this file, so it includes a lot of different
    stuff, some of which was never used for anything important.  However, many functions here are useful.
  
  Aij_validation_driver.R : Driver script for Aij_validation.R
  
    I almost didn't include this, but this is my (mess of a) driver script for generating a lot of the 
    statistics and graphs from Aij_validation.R at once.  Stuff is commented in or out as necessary.
    Basically this code is a gigantic example of how to use Aij_validation.R. 
  
  outliers.R : Driver script for examining differences in outlier handling
  
  demo_simulation.R : This script deals with examining the operon-confidence enhancement for Bayesian 2.0.
    (see the relevant submitted paper).  It uses a different data-simulation methodology. 
  
-------------------------------------------------------
The following R files relate to the simulation of data for the bacterial project.

  simulate_data.R : Contains the primary function for generating simulated data, based on
    many different parameters
    
  simulation_from_fits.R : Contains the code required for generating simulated data based 
    on the Gibbs-sampler-fitted distributions drawn from some real expression data

-------------------------------------------------------
The following are other miscellaneous R files related to the Bayesian project.

  bayesfactorplots.R : Functions I used to generate various graphs for illustrative purposes,
    and investigate the scaling behavior of various Bayes-Factor-y functions on various kinds
    of simulated data. Apologies for the poor documentation. "LBF" means "Ln(Bayes Factor)" and
    "NBF" means "Normalized Bayes Factor". Note this NBF is the univariate or "old method" NBF
    which doesn't correct for operon size. 
  
  Clean_matt_calls_April_2016.R : When Matt DeJongh inevitably sends newer gold-standard calls based
    on FVA, this is the script to use to standardize the data format to what our scripts expect.  
    At least, this worked for the calls he sent in April 2016. 
    
  ishii_analysis_script.R : File for analysis of the "Ishii data" (see the Bayesian 1.0 publication).
    This borrows heavily from the iMRM code (see below) because of the need to apply GPR relationships. 

-------------------------------------------------------
The following R files relate to Aim 2 specifically. (I think.)

  pilot_aim2.R : Functions relating to the Aim 2 "pilot study" (which isn't really pilot anymore,
    but the file didn't get renamed) on simulated data.  This script was used in generating the 
    "Aim 2 new analysis 2-3-16.docx" document, which explains things pretty thoroughly. 
  
  tempaim2script.R : Contains a variety of functions related to data analysis (on real data)
    for Aim 2.  Some of these are more like utility functions for Aim 2-specific uses.
    (Originally a temporary script which gradually morphed into the code you see now, and was
    never renamed.)
  
  analyze_scenarios_complete.R : Functions relating to scenarios (pathways)
  
  nonCompleteExpts_hope_ecoli.Rdata : Data file used by analyze_scenarios_complete (and maybe
    some other scripts) which is simply a list (vector technically) of the experiments in the
    Hope Ecoli dataset which are not complete, that is, do not contain yeast extract.
  
  explore_associations.R : Systematic analysis of the media present in each experiment
    including significance tests (Fisher's exact test) for the association between various 
    "features" and gene activity states. 
  
-------------------------------------------------------
The following R files relate to the iMRM project specifically.

  iMRM.R : The main iMRM file, which contains the core functions implementing the iMRM method.
    Calls functions found in many of the other files. 
  
  iMRM CPLEX.R : Code for interfacing with the CPLEX solver (and setting up the MILP itself). 
  
  iMRM gene expression and kos.R : Code for computing the penalties based on gene expression,
    and handling gene knockouts.  Called by iMRM CPLEX.R.
    
  iMRM media.R : Code for dealing with varied media conditions, by a couple different methods.
  
  iMRM driver.R : Driver script for the iMRM code; also contains all of the code for parsing and 
    loading input from our input data files. 
  
  iMRM gene calls.R : Now-abandoned code attempting to produce "gold standard" calls based on 
    running the iMRM with no gene expression input.  We basically decided the FVA method (which
    Matt DeJongh typically runs for us) is superior.  Nonetheless, I'm including the code in 
    case this ever gets revisited. 

-------------------------------------------------------
----------*** DOCUMENTATION AND WRITE-UPS ***----------
-------------------------------------------------------


-------------------------------------------------------
Bayesian 1.0 and 2.0 related

Primarily I encourage you to read the submitted conference paper for Bayesian 2.0, and as 
necessary, the published paper for Bayesian 1.0. (Ask Dr. Tintle for final copies of these.)

I have also included a very early draft of the Bayesian 2.0 paper which included many more 
details which were cut out later for length reasons. (And lots of comments, some of which
might serve as a basis for inspiration for improvements in Bayesian 3.0 if you are looking
for ideas.)  

-------------------------------------------------------
Aim 2 related

I've included the latest paper outline (from June 2016) for Aim 2, which describes the
current vision for an in silico Aim 2 paper.  

Also the latest "results" on Carrera and M3D data, together with "Explanation of scores.docx" 
which provides some interpretation of the data there. Note these may now be outdated and not 
based on the final NBF-based Bayesian 2.0 aijs; you can regenerate this yourself if you want 
with some of the functions in tempaim2script.R. 

-------------------------------------------------------
iMRM related

Most important here is the "iMRM method description.docx" - latest version January 2016. This
will get you pretty much up to speed on the iMRM project and explains almost all of the code.

The "thoughts on maxPenalty" are some musings of mine which you may or may not find useful. 

I've also included a nifty little graphic I made in August 2015; it compares the current iMRM
procedure (pretty much still describes the current procedure, there haven't been major changes)
with a competing procedure outlined in the Carrera paper (find it online in the journal Molecular
Systems Biology with DOI 10.15252/msb.20145108, published May 2014).  It also proposes a "merger"
of the two methods which takes (in my opinion) the best from both methods.  Specifically it would
incorporate more data sources than the current iMRM, and also remove the restriction that the iMRM
can only make reaction-flux predictions for experiments which have actually physically been performed,
as opposed to being able to make in silico predictions about hypothetical experiments, as the 
Carrera method can.  No code has been written for this integration, it's more a crazy brainchild of mine.

-------------------------------------------------------
---------------------*** DATA ***----------------------
-------------------------------------------------------


I've included all our current input files in the inputs/ folder.  No other explanation given here;
if you're wondering about a particular file, just find where it's mentioned in the code (Linux's
grep tool is awesome).  

I've included some cache files in the Rsave/ directory.  For space reasons I didn't include everything.
This means you'll have to spend a while regnerating files for whatever you decide you want.  You should
probably run things on Beaker when doing anything computationally intensive.  That said, if there are a
couple more specific files you want, just let me know and I'll post them to the shared folder. 