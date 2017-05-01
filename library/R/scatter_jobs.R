write_R_files = function( bsub_argv, var_bsub_argv, control_argv )
{
   R_DIR = control_argv$R_DIR

   cat("Creating R files")
#   for( i in 1:nrow(var_bsub_argv) ) 
   res = foreach( i = 1:nrow(var_bsub_argv) ) %dopar%
   {
      for( j in 1:ncol(var_bsub_argv) )
      {
         bsub_argv[[ colnames(var_bsub_argv)[j] ]] = var_bsub_argv[i,j] 
      }

      TMPFN = tempfile(pattern="tmp", tmpdir="." )
      for( j in 1:length(bsub_argv) )
      {
         A = ifelse( j == 1, FALSE, TRUE ) ## overwrite exiting file
         if( is.character( bsub_argv[[j]] ) )
         {
            cat( names(bsub_argv)[j], "= \"",  bsub_argv[[j]], "\"\n", file=TMPFN, sep="", append=A)
         }
         else
         {
            cat( names(bsub_argv)[j], "= ",  bsub_argv[[j]], "\n", file=TMPFN, sep="", append=A )
         }
      } 
      
     ## create .R executable code
      R_FN = paste( i, ".R", sep="")
      qjr = file.path( R_DIR, R_FN )
      sc = paste( "cat ", TMPFN, " ", control_argv$R_STUB_FN, " > ", qjr, sep="" )
      system(sc)
      
      sc = paste( "rm ", TMPFN, sep="")
      system(sc)
      cat(".")
      return()
   }
   cat("done\n")
}

# FIXME
# Dispatch now works, but we are not collecting the job ids,
# which causes failure in scatter_jobs
# Since multiple job ids will need to be monitored, the current code will
# *not* for SGE.
# Task array dispatch seems to work on SGE (invoked by specifycing EXE_ENGINE
# as UGER)... perhaps UGER and SGE should be consolidated?
# Should serial dispatch still be supported?
serial_dispatch = function( control_argv, var_bsub_argv )
{
	R_DIR = control_argv$R_DIR
# now dispatch jobs
   res = foreach( i = 1:nrow(var_bsub_argv) ) %dopar%
   {
       R_FN = paste( i, ".R", sep="")
       jname = paste( R_DIR, i, sep="" )
       qjname = paste( "\"", jname, "\"", sep="" )
       qjout = paste( "\"", jname, ".bsub.out", "\"", sep="" )
       qjerr = paste( "\"", jname, ".bsub.err", "\"", sep="" )
       qjrout = paste( jname, ".R.out", sep="" )
       qjr = paste( R_DIR, R_FN, sep="" )

       bsub = paste( "R < ", qjr, " --no-restore --no-save > ", qjrout , sep="" )
#      cat( bsub, file="bsub.tmp", append=F )

#      prescript = paste( "\"ls ", BASE_DIR, "\"", sep="" )
#      sc = paste( "bsub < bsub.tmp -r -mig 5 -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", bjob, " -E ", prescript, sep="" )

      if( control_argv[["EXE_ENGINE"]] == "LSF" )
      {
         if( control_argv$QUEUE %in% c("hour", "bhour") )
         {
#            sc = paste( "bsub < bsub.tmp -W 4:00 -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", control_argv$BJOB, sep="" )
            sc = paste(control_argv$engine_setup, " echo \"", bsub, "\" | bsub -W 4:00 -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", control_argv$BJOB, sep="" )
         } 
         else {
#            sc = paste( "bsub < bsub.tmp -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", control_argv$BJOB, sep="" )
            sc = paste(control_argv$engine_setup, " echo \"", bsub, "\" | bsub -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", control_argv$BJOB, sep="" )
         } 
         if( control_argv[["groupname"]] != "" )
         {
            sc = paste(sc, " -G ", control_argv[["groupname"]], sep="")
         }

      }

      if( control_argv[["EXE_ENGINE"]] == "SGE" )
      {
#            sc = paste( "bsub < bsub.tmp -q ", control_argv$QUEUE, " -o ", qjout, " -e ",  qjerr, " -J ", control_argv$BJOB, sep="" )
         sc = paste(control_argv$engine_setup,  " echo \"", bsub, "\" | qsub -o ", qjout, " -e ", qjerr, " -N ", control_argv$BJOB, sep="" )
       
         if( control_argv[["groupname"]] != "" )
         {
            sc = paste(sc, " -G ", control_argv[["groupname"]], sep="")
         }
         if( control_argv[["QUEUE"]] != "" )
         {
            sc = paste( sc, " -q ", control_argv[["QUEUE"]], sep="")
         }
      }

      if( control_argv[["EXE_ENGINE"]] == "MULTICORE" )
      { sc = bsub }

      if( !control_argv$DRY_RUN ) 
      {
         system(sc)        
      }
      return()
   }

   if( control_argv$wait==TRUE )
   {
      if( control_argv$EXE_ENGINE %in% c("LSF") )
      {
         while ( system(paste("bjobs -J ", control_argv$BJOB, " | wc -l"), intern=T)>0 )
         { Sys.sleep(60) }
      }
   }
}


scatter_jobs = function( control_argv, bsub_argv, var_bsub_argv )
{
   if( !(control_argv[["EXE_ENGINE"]] %in% c("LSF", "SGE", "UGER", "MULTICORE")) ) {stop("Unsupported EXE_ENGINE specified in control_argv: must be one of (LSF, SGE, UGER, MULTICORE)") }

   R_DIR = control_argv$R_DIR
   dir.create( R_DIR, recursive = TRUE, showWarnings = FALSE)

   write_R_files(bsub_argv, var_bsub_argv, control_argv)

   if( !control_argv$DRY_RUN & control_argv[["EXE_ENGINE"]] %in% c("LSF", "SGE", "MULTICORE")) 
   { 
      serial_dispatch(control_argv, var_bsub_argv )
   }


## for UGER, submit all jobs as one using a task array
   if( !control_argv$DRY_RUN & control_argv[["EXE_ENGINE"]] == "UGER" ) 
   {
      N_tasks = nrow(var_bsub_argv)
      R_FN = paste( "$SGE_TASK_ID", ".R", sep="")
      jname = paste( R_DIR, "$SGE_TASK_ID", sep="" )
      qjname = paste( "\"", jname, "\"", sep="" )
 #     qjout = paste( "\"", jname, ".bsub.out", "\"", sep="" )
#      qjerr = paste( "\"", jname, ".bsub.err", "\"", sep="" )
      qjrout = paste( jname, ".R.out", sep="" )
      qjr = paste( R_DIR, R_FN, sep="" )

# >& redirects both stdout and stderr
      bsub = paste( "R < ", qjr, " --no-restore --no-save >& ", qjrout , sep="" )

      SH_FN = paste( control_argv$BJOB, ".qsub.sh", sep="" )
      cat("#!/bin/bash\n", file=SH_FN, append=FALSE)
      cat(bsub, file=SH_FN, append=TRUE)

      sc = paste(control_argv$engine_setup, "qsub", "-q", control_argv[["QUEUE"]], sprintf("-t 1-%d", N_tasks), "-tc", N_tasks, "-cwd", "-V", "-o /dev/null", "-e /dev/null", "-l h_vmem=4g", "-N", control_argv$BJOB, SH_FN)
      
      stdout = system(sc, intern=TRUE)
# extract job id from qsub stdout
      #job.id = strsplit(stdout, "\\." )
      #job.id = strsplit(job.id[[2]], "Your job-array " )
      #job.id = job.id[[1]][2]
			stdout = paste(stdout, collapse="\n")
      job.id = sub(".*Your job-array ([0-9]+).*", "\\1", stdout)
      print( paste("Captured qsub job.id ", job.id, sep=""))

      if(is.na(job.id)) { stop(paste("Invalid job.id.  Captured: ", stdout, " from stdout.", sep="") ) }
   }


   if( !control_argv$DRY_RUN & control_argv$wait==TRUE )
   {
      if( control_argv$EXE_ENGINE %in% c("UGER","SGE") )
      {
         while ( TRUE )
         {
            stdout = system(paste(control_argv$engine_setup, "qstat"), intern=TRUE )
            # FIXME if EXE_ENGINE == "SGE", job.id will be undefined here!
            if( length( grep( job.id, stdout, value=TRUE ) ) == 0 ) { break }
            Sys.sleep(60) 
         }

       ## if the short queue was used, look for tasks that were killed for going over the 2hr (7200 second) time limit
         if( control_argv[["QUEUE"]] == "short" ) 
         {
            user.id = Sys.getenv("USER")
            cmd = paste( "cat /broad/uge/research/research/common/accounting | grep ", user.id, " | grep short | grep ", job.id, " | awk -F: -v max=7198 \' $14 >= max {print $36}\'", sep="")
# {print $6\".\"$36\" : \"$14}\'", sep="")
            stdout = system(cmd, intern=TRUE)
            failed = as.integer(stdout)
            if(length(failed) == 0 ) { return() }

            print( paste("Detected ", length(failed), " of ", N_tasks, " exceeded 2hr run limit on the short queue; resubmitting to the long queue.", sep=""))
       ## now resubmit them to the long queue
            N_failed = length(failed)
            var_bsub_argv = var_bsub_argv[ failed,, drop=FALSE ]

           # re-write R-files for failed jobs into a 'long' subdir
            R_DIR = file.path(control_argv$R_DIR, "long")
            dir.create( R_DIR, recursive = TRUE, showWarnings = FALSE)
            control_argv[["R_DIR"]] = R_DIR
            write_R_files(bsub_argv, var_bsub_argv, control_argv)

          ## repoint exec script to long subdir
            qjr = file.path( R_DIR, R_FN )
            bsub = paste( "R < ", qjr, " --no-restore --no-save >& ", qjrout , sep="" )
            SH_FN = paste( control_argv$BJOB, ".qsub.sh", sep="" )
            cat("#!/bin/bash\n", file=SH_FN, append=FALSE)
            cat(bsub, file=SH_FN, append=TRUE)

            sc = paste(control_argv$engine_setup, "qsub", "-q long", sprintf("-t 1-%d", N_failed), "-tc", N_failed, "-cwd", "-V", "-o /dev/null", "-e /dev/null", "-N", control_argv$BJOB, SH_FN)
      
            stdout = system(sc, intern=TRUE)
            # extract job id from qsub stdout
            #job.id = strsplit(stdout, "\\." )
            #job.id = strsplit(job.id[[2]], "Your job-array " )
            #job.id = job.id[[1]][2]
						stdout = paste(stdout, collapse="\n")
            job.id = sub(".*Your job-array ([0-9]+).*", "\\1", stdout)
            print( paste("Captured qsub job.id ", job.id, sep=""))

            if(is.na(job.id)) { stop(paste("Invalid job.id.  Captured: ", stdout, " from stdout.", sep="") ) }

            while ( TRUE )
            {
               stdout = system(paste(control_argv$engine_setup, "qstat"), intern=TRUE )
               if( length( grep( job.id, stdout, value=TRUE ) ) == 0 ) { break }
               Sys.sleep(60) 
            }
         }
      }
   }
}


