#Optional input function. Counts sequences and creates K count-sequence pairs. Provides input variable to thresholdSeq().
#Arguments: 
#a -- A vector of sequences/features.
#Returns: 
#(unnamed) -- A vector of count values.
makeCountSequencePairs <- function(a) {
    return(as.numeric(table(a))); 
}

#Randomly resample the K unique count-sequence pairs. Used in the first step of the Threshold-seq algorithm. 
#Arguments:
#d -- the original set of K count-sequence pairs
#n -- the number of items to resample: SHOULD ALWAYS BE EQUAL TO K
#Returns:
#(unnamed) -- A vector of K resampled count values. 
randPerm <- function(d,n) {
    return(sample(d,n,replace=TRUE));	
}

#Find the mode of a numerical vector. Used to determine the final threshold of the Threshold-seq algorithm. 
#Arguments:
#a -- A numerical vector. 
#Returns:
#(unnamed) -- An integer as the mode of the input argument. 
nmode <- function(a) {
    ux <- unique(a)
    return(ux[which.max(tabulate(match(a, ux)))])	
}

#Find the threshold for this dataset.
#Arguments:
#b -- table containing x values for each of the N iterations at each CDF value from MIN to MAX
#c -- all possible alpha values from MIN to MAX, stepping by step
#Returns:
#out -- An object containing the CDFtarget value (CDFtarget) and the threshold for this dataset (threshold). 
findThresh <- function(b,c) {
    
    out=NULL;
    
    #determine the CDFtarget value and threshold
    for (i in 1:nrow(b)) {
        if (length(unique(b[i,]))>1 && max(table(b[i,]))/sum(table(b[i,])) < 0.99) {
            out$CDFtarget=c[i]; 
            out$threshold=nmode(b[i,]);
            break;
        }	
    }
    
    #return the output object
    return(out);
}

#Find the read threshold for this sample.
#Arguments:
#d -- The input data of K count-sequence pairs. Here, a numerical vector containing the counts for the K unique sequences.
#output -- The text file to which to write the threshold output value. 
#nperm -- The number N of iterations through the Threshold-seq algorithm to run.
#CDFmin -- The minimum value of the CDF to begin the scan from. 
#CDFstep -- The value by which to increase the scanned CDF values. 
#verbose -- Should incremental output be typed to the console?
#autoRetry -- If a threshold can't be found, should Threshold-seq automatically reattempt with a CDFstep/2? 
#Returns:
#out -- An object containing the threshold for this dataset (threshold), the table of x values determined over the range of CDF values scanned (thresholdTable), the CDFtarget value (CDFtarget), the original data (d), and the original data CDF (d_CDF). 
thresholdSeq <- function(d,nperm=1000,CDFmin=0.90,CDFstep=0.005,verbose=TRUE,autoRetry=TRUE) {
    #Bring in the data in the d object. This should be a numerical vector containing the count values corresponding to K count-sequence pairs. 
    if (length(d) == 0) { print('Error in input. Exiting. Please supply at least: \n d - vector of integers representing read counts to be thresholded \n and output - filename to which threshold should be appended.'); return(NULL); }	
    if (verbose) { print('Input data summary:--'); print(summary(d)); print('---'); }
    
    sanityCheck = as.numeric(d);
    if (length(which(is.na(sanityCheck)))>0) {
        print('Error: non-numeric data found in input. Please double check your input and relaunch Threshold-seq.');
        write('Error: non-numeric data found in input. Please double check your input and relaunch Threshold-seq.',file=output,append=TRUE);
        return(NULL);
    }
    
    #Establish the CDF values to be scanned from CDFmin to 1.0, stepping by CDFstep.
    alpha=seq(CDFmin,0.99,CDFstep);
    
    #Build the necessary storage variables. 
    thresholds=matrix(0,length(alpha),nperm);
    
    #Iterate over all CDF values and permutations 
    for (s2 in 1:nperm) {
        #First, randomly resample the data.
        CDF=ecdf(randPerm(d,length(d)));
        #precompute all the CDF values 
        unq_d=sort(unique(d));
        CDFvals=CDF(unq_d);
        
        for (s1 in 1:length(alpha)) {
            ind=which(CDFvals>=alpha[s1])[1];
            thresholds[s1,s2]=unq_d[ind];  
        }
    }
    
    #Find the threshold based on the information gathered above. 
    threshInfo=findThresh(thresholds,alpha);
    
    #If we weren't able to find a threshold with these features, warn the user    and start over if the flag is set. 
    if (length(threshInfo$thresh) == 0 && autoRetry == TRUE) {
        write(paste('Unable to find a threshold with parameters:\nnperm=',nperm,'\nCDFmin=',CDFmin,'\nCDFstep=',CDFstep,'\nReattempting with new CDFstep=',CDFstep/2,sep=''),file=output,append=TRUE);
        return(thresholdSeq(d,output,nperm,CDFmin,CDFstep/2,verbose,autoRetry));
    }
    
    #Print threshold if verbosity set. 
    if (verbose) { 	
        print(paste('Read Threshold=',threshInfo$thresh)); 
        print(paste('CDFtarget=',threshInfo$CDFtarget));
    }
    
    #Find and return the original CDF of the data. 
    CDFd=ecdf(d);
    d_CDF=cbind(unique(d),CDFd(unique(d)));
    d_CDF=d_CDF[order(d_CDF[,1]),];
    
    #Write the threshold as a single number to an output file.
    # write(paste('Threshold-seq threshold = ',as.character(threshInfo$thresh),sep='\t'),file=output,append=TRUE);
    
    #Return the output object. 
    out=NULL;
    out$threshold=threshInfo$thresh; 
    out$thresholdTable=thresholds;
    out$CDFtarget=threshInfo$CDFtarget; 
    out$d=d;
    out$d_CDF=d_CDF;
    
    return(out); 
}
