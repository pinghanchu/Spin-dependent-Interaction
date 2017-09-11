#!/usr/bin/perl
#my $scriptpath = "\$GATDIR/MJDCalibration/";
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/Spin-dependent-Interaction/";
my $cal = $scriptpath."spindependentinteraction";
my $plot = $scriptpath."plot";
my @force =    (2,3,4,6,8,9,11,12,14,15,16);
my @position = (2,2,2,1,1,2, 2, 1, 2, 2, 2);
my @max = (-10,0,0,-2,-2,-10,-10,-10,-10,0,3);
my @min =(-40,-30,-40,-18,-18,-40,-30,-35,-30,-20,-7);
my $size = @force;
for(my $i = 0;$i<$size;$i++){
    my $iforce = $force[$i];
    my $ipos = $position[$i];
    my $imax = $max[$i];
    my $imin = $min[$i];
    my $file = $iforce."_".$ipos;
    my $app  = "cal_".$file.".csh";
    print "submit ", $app,"\n";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "$cal $iforce $ipos\n";
    print $fh "$plot $iforce $ipos $imax $imin\n";
    close $fh;
    system("chmod 755 $app");
    system("./$app");
    #system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
    #system("$calibration $startrun $endrun $coverstartrun $coverendrun $ienr $ipos");

}







