#!/usr/bin/perl -w

## This perl script calls feidtrs_joao.exe for all of the .csv files 
#in the current directory and saves the results in .csv files by prefixinf
#the original name by prefixing d_indices_. The indices used are specified 
#by the list below

use File::Find;

$crystal = 6; 
@indices= ('1 0 -1 0', '0 0 0 2', '1 0 -1 1', '1 0 -1 2', '1 1 -2 0', '1 0 -1 3', '1 1 -2 2', '2 0 -2 1' ); 
$ca=1.624; #Mg
#$ca=1.587; #Ti
#$ca=1.593; #Zr

$angle_tol=10;

foreach (@indices) {
    $nice_indices=$_;
    $nice_indices=~s/\s//g;
    $indices=$_;
    print "\nReflection: $_ \n";

    find(\&wanted,'.');


    sub wanted { 
        if (m/^[^d]\S+\.csv/) {
        open (TMP_FILE,">","tmp") ||die "can't open tmp file $!"; 
        $file_in=$_;
        $file_out="d\_$nice_indices\_$_";
        print TMP_FILE "$crystal\n$ca\n$indices\n$angle_tol\n$file_in\n$file_out\n";
        close(TMP_FILE) || die "can't close tmp file"; 
        #print "Doing $_:";
	 system("/home/mclssgt2/bin/feistrs < tmp > log");
        #print " done\n";

        }
    }
}   
  
