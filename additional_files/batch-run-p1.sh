rm *.csv *.dat *.frs;

ls -latr /home/mclssgt2/bin > verinfo
/home/mclssgt2/bin/fs3pre < fs3pre_inputs.txt ;
/home/mclssgt2/bin/fasolt3 < fe_input2.txt ;
# perl /home/tsivoulas/bin/perl/gstress.pl <fafrep.txt;
# perl /home/tsivoulas/bin/perl/csvmkr_1.pl;
# perl /home/tsivoulas/bin/perl/feistrs_do_fcc.pl;
# perl /home/tsivoulas/bin/perl/extract_means_fcc.pl;
