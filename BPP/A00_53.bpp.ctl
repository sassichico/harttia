    seed =  54321

    #seqfile = ./BPP/Harttia/Subsets/Harttia_31_bpp.txt
    #Imapfile = ./BPP/Harttia/Subsets/Harttia_Imap.txt
    #outfile = ./BPP/Harttia/Output31/Harttia_31_out.txt
    #mcmcfile = ./BPP/Harttia/Output31/Harttia_31_mcmc.txt

    seqfile = /home/fernando/BPP/Harttia/Subsets/Harttia_51_bpp.txt
    Imapfile = /home/fernando/BPP/Harttia/Subsets/Harttia_Imap.txt
    outfile = ./Harttia_53_out.txt
    mcmcfile = ./Harttia_53_mcmc.txt

#Optional
    constraintfile = /home/fernando/BPP/Harttia/Subsets/constraints.txt

    # fixed species tree and delimitation
    speciesdelimitation = 0 
    speciestree = 0

    species&tree = 18 villasboas rondoni dissidens duruventris spRP punctata carvalhoi spRM torrenticola intermontana longipinna gracilis loricariformis kronei guianensis fowleri spBG farlowella 
                   6 6 6 6 6 6 6 5 6 6 6 6 6 6 6 6 6 6
                   (((((kronei,(loricariformis,(gracilis,(longipinna,(intermontana,(torrenticola,(spRM,carvalhoi))))))),(((duruventris,((villasboas,rondoni),dissidens)),spRP),punctata)),(guianensis,fowleri)),spBG),farlowella);
    # sequence data are unphased for all 4 populations
    phase =   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1

    model = HKY
    # 0: no data (prior); 1:seq like
    usedata = 1  

    # number of data sets in seqfile
    nloci =   500
    #nloci =   179

    # remove sites with ambiguity data (1:yes, 0:no)?
    cleandata = 0    

    # invgamma(a, b) for theta
    thetaprior = invgamma 3 0.002  

    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's
    tauprior = invgamma 3 0.18    

    *     heredity = 0
    *     locusrate = 0

    # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
    finetune =  1

    # MCMC samples, locusrate, heredityscalars, Genetrees
    #print = 1 1 1 1 1 
    #print = -1 0 0 0 0    
    print = 1 0 0 0 1

    #Checkpoint files
    checkpoint = 8000 8000

    #The total number of MCMC iterations is burnin + nsample × sampfreq
    burnin = 40000

    sampfreq = 2

    #nsample = 50000
    nsample = 400000

    threads = 22
