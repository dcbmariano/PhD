#!/usr/bin/perl -w 
# *******************************************************
# *		Structural Bioinformatics Group		*
# *		UFMG/DCC				*
# *   ----------------------------------------------	*
# *							*
# * Douglas Eduardo Valente Pires - dpires@dcc.ufmg.br	*
# * www.dcc.ufmg.br/~dpires				*
# *   ----------------------------------------------	*
# *							*
# * Modificado por Diego Mariano                        *
# * Posicao das coordenadas x, y e z foram modificadas  *
# *******************************************************

use strict;
use warnings;

sub trim;
sub distance;

# ____________________________________________________________________________________________________________________
my %Hydrophobic;
$Hydrophobic{"ALACB"} = 1;
$Hydrophobic{"ARGCB"} = $Hydrophobic{"ARGCG"} = $Hydrophobic{"ARGCD"} = 1;
$Hydrophobic{"ASNCB"} = 1;
$Hydrophobic{"ASPCB"} = 1;
$Hydrophobic{"CYSCB"} = 1;
$Hydrophobic{"GLNCB"} = $Hydrophobic{"GLNCG"} = 1;
$Hydrophobic{"GLUCB"} = $Hydrophobic{"GLUCG"} = 1;
$Hydrophobic{"HISCB"} = $Hydrophobic{"HISCG"} = $Hydrophobic{"HISCD2"} = $Hydrophobic{"HISCE1"} = 1;
$Hydrophobic{"ILECB"} = $Hydrophobic{"ILECG1"} = $Hydrophobic{"ILECG2"} = $Hydrophobic{"ILECD1"} = 1;
$Hydrophobic{"LEUCB"} = $Hydrophobic{"LEUCG"} = $Hydrophobic{"LEUCD1"} = $Hydrophobic{"LEUCD2"} = 1;
$Hydrophobic{"LYSCB"} = $Hydrophobic{"LYSCG"} = $Hydrophobic{"LYSCD"} = 1;
$Hydrophobic{"METCB"} = $Hydrophobic{"METCG"} = $Hydrophobic{"METCE"} = 1;
$Hydrophobic{"PHECB"} = $Hydrophobic{"PHECG"} = $Hydrophobic{"PHECD1"} = $Hydrophobic{"PHECD2"} = $Hydrophobic{"PHECE1"} = $Hydrophobic{"PHECE2"} = $Hydrophobic{"PHECZ"} = 1;
$Hydrophobic{"PROCB"} = $Hydrophobic{"PROCG"} = $Hydrophobic{"PROCD"} = 1;
$Hydrophobic{"THRCG2"} = 1;
$Hydrophobic{"TRPCB"} = $Hydrophobic{"TRPCG"} = $Hydrophobic{"TRPCD1"} = $Hydrophobic{"TRPCD2"} = $Hydrophobic{"TRPCE2"} = $Hydrophobic{"TRPCE3"} = $Hydrophobic{"TRPCH2"} = $Hydrophobic{"TRPCZ"} = $Hydrophobic{"TRPCZ2"} = $Hydrophobic{"TRPCZ3"} = 1;
$Hydrophobic{"TYRCB"} = $Hydrophobic{"TYRCG"} = $Hydrophobic{"TYRCD1"} = $Hydrophobic{"TYRCD2"} = $Hydrophobic{"TYRCE1"} = $Hydrophobic{"TYRCE2"} = $Hydrophobic{"TYRCZ"} = 1;
$Hydrophobic{"VALCB"} = $Hydrophobic{"VALCG1"} = $Hydrophobic{"VALCG2"} = 1;

my %Positive;
$Positive{"ARGNH1"} = $Positive{"ARGNH2"} = 1;
$Positive{"HISND1"} = $Positive{"HISNE2"} = 1;
$Positive{"LYSNZ"} = 1;

my %Negative;
$Negative{"ASPOD1"} = $Negative{"ASPOD2"} = 1;
$Negative{"GLUOE1"} = $Negative{"GLUOE2"} = 1;

my %Acceptor;
$Acceptor{"ALAO"} = 1;
$Acceptor{"ARGO"} = 1;
$Acceptor{"ASNO"} = $Acceptor{"ASNOD1"} = 1;
$Acceptor{"ASPO"} = $Acceptor{"ASPOD1"} = $Acceptor{"ASPOD2"} = 1;
$Acceptor{"CYSO"} = 1;
$Acceptor{"GLNO"} = $Acceptor{"GLNOE1"} = 1;
$Acceptor{"GLUO"} = $Acceptor{"GLUOE1"} = $Acceptor{"GLUOE2"} = 1;
$Acceptor{"GLYO"} = 1;
$Acceptor{"HISO"} = 1;
$Acceptor{"ILEO"} = 1;
$Acceptor{"LEUO"} = 1;
$Acceptor{"LYSO"} = 1;
$Acceptor{"METO"} = 1;
$Acceptor{"PHEO"} = 1;
$Acceptor{"PROO"} = 1;
$Acceptor{"SERO"} = 1;
$Acceptor{"THRO"} = 1;
$Acceptor{"TRPO"} = 1;
$Acceptor{"TYRO"} = 1;
$Acceptor{"VALO"} = 1;

my %Donor;
$Donor{"ALAN"} = 1;
$Donor{"ARGN"} = $Donor{"ARGNE"} = $Donor{"ARGNH1"} = $Donor{"ARGNH2"} = 1;
$Donor{"ASNN"} = $Donor{"ASNND2"} = $Donor{"ASNOD1"} = 1;
$Donor{"ASPN"} = 1;
$Donor{"CYSN"} = 1;
$Donor{"GLNN"} = $Donor{"GLNNE2"} = 1;
$Donor{"GLUN"} = 1;
$Donor{"GLYN"} = 1;
$Donor{"HISN"} = $Donor{"HISND1"} = $Donor{"HISNE2"} = 1;
$Donor{"ILEN"} = 1;
$Donor{"LEUN"} = 1;
$Donor{"LYSN"} = $Donor{"LYSNZ"} = 1;
$Donor{"METN"} = 1;
$Donor{"PHEN"} = 1;
$Donor{"PRON"} = 1;
$Donor{"SERN"} = $Donor{"SEROG"} = 1;
$Donor{"THRN"} = $Donor{"THROG1"} = 1;
$Donor{"TRPN"} = $Donor{"TRPNE1"} = 1;
$Donor{"TYRN"} = $Donor{"TYROH"} = 1;
$Donor{"VALN"} = 1; 

my %Aromatic;
$Aromatic{"HISCG"} = $Aromatic{"HISND1"} = $Aromatic{"HISCD2"} = $Aromatic{"HISCE1"} = $Aromatic{"HISNE2"} = 1;
$Aromatic{"PHECG"} = $Aromatic{"PHECD1"} = $Aromatic{"PHECD2"} = $Aromatic{"PHECE1"} = $Aromatic{"PHECE2"} = $Aromatic{"PHECZ"} = 1;
$Aromatic{"TRPCG"} = $Aromatic{"TRPCD1"} = $Aromatic{"TRPCD2"} = $Aromatic{"TRPNE1"} = $Aromatic{"TRPCE2"} = $Aromatic{"TRPCE3"} = $Aromatic{"TRPCZ2"} = $Aromatic{"TRPCZ3"} = $Aromatic{"TRPCH2"} = 1;
$Aromatic{"TYRCD1"} = $Aromatic{"TYRCD2"} = $Aromatic{"TYRCE1"} = $Aromatic{"TYRCE2"} = $Aromatic{"TYRCG"} = $Aromatic{"TYRCZ"} = 1;

my %Sulfur;
$Sulfur{"CYSS"} = 1; $Sulfur{"METSD"} = 1;

# ____________________________________________________________________________________________________________________
my %AtomFeature;
$AtomFeature{"ALACB"} =  "Hydro";
$AtomFeature{"ARGCB"} = "Hydro";
$AtomFeature{"ARGCG"} = "Hydro";
$AtomFeature{"ARGCD"} =  "Hydro";
$AtomFeature{"ASNCB"} =  "Hydro";
$AtomFeature{"ASPCB"} =  "Hydro";
$AtomFeature{"CYSCB"} =  "Hydro";
$AtomFeature{"GLNCB"} = "Hydro";
$AtomFeature{"GLNCG"} =  "Hydro";
$AtomFeature{"GLUCB"} = "Hydro";
$AtomFeature{"GLUCG"} =  "Hydro";
$AtomFeature{"HISCB"} = "Hydro";
$AtomFeature{"HISCG"} = "Hydro";
$AtomFeature{"HISCD2"} = "Hydro";
$AtomFeature{"HISCE1"} =  "Hydro";
$AtomFeature{"ILECB"} = "Hydro";
$AtomFeature{"ILECG1"} = "Hydro";
$AtomFeature{"ILECG2"} = "Hydro";
$AtomFeature{"ILECD1"} =  "Hydro";
$AtomFeature{"LEUCB"} = "Hydro";
$AtomFeature{"LEUCG"} = "Hydro";
$AtomFeature{"LEUCD1"} = "Hydro";
$AtomFeature{"LEUCD2"} =  "Hydro";
$AtomFeature{"LYSCB"} = "Hydro";
$AtomFeature{"LYSCG"} = "Hydro";
$AtomFeature{"LYSCD"} =  "Hydro";
$AtomFeature{"METCB"} = "Hydro";
$AtomFeature{"METCG"} = "Hydro";
$AtomFeature{"METCE"} =  "Hydro";
$AtomFeature{"PHECB"} = "Hydro";
$AtomFeature{"PHECG"} = "Hydro";
$AtomFeature{"PHECD1"} = "Hydro";
$AtomFeature{"PHECD2"} = "Hydro";
$AtomFeature{"PHECE1"} = "Hydro";
$AtomFeature{"PHECE2"} = "Hydro";
$AtomFeature{"PHECZ"} =  "Hydro";
$AtomFeature{"PROCB"} = "Hydro";
$AtomFeature{"PROCG"} = "Hydro";
$AtomFeature{"PROCD"} =  "Hydro";
$AtomFeature{"THRCG2"} =  "Hydro";
$AtomFeature{"TRPCB"} = "Hydro";
$AtomFeature{"TRPCG"} = "Hydro";
$AtomFeature{"TRPCD1"} = "Hydro";
$AtomFeature{"TRPCD2"} = "Hydro";
$AtomFeature{"TRPCE2"} = "Hydro";
$AtomFeature{"TRPCE3"} = "Hydro";
$AtomFeature{"TRPCH2"} = "Hydro";
$AtomFeature{"TRPCZ"} = "Hydro";
$AtomFeature{"TRPCZ2"} = "Hydro";
$AtomFeature{"TRPCZ3"} =  "Hydro";
$AtomFeature{"TYRCB"} = "Hydro";
$AtomFeature{"TYRCG"} = "Hydro";
$AtomFeature{"TYRCD1"} = "Hydro";
$AtomFeature{"TYRCD2"} = "Hydro";
$AtomFeature{"TYRCE1"} = "Hydro";
$AtomFeature{"TYRCE2"} = "Hydro";
$AtomFeature{"TYRCZ"} =  "Hydro";
$AtomFeature{"VALCB"} = "Hydro";
$AtomFeature{"VALCG1"} = "Hydro";
$AtomFeature{"VALCG2"} =  "Hydro";

$AtomFeature{"ARGNH1"} = "Pos";
$AtomFeature{"ARGNH2"} = "Pos";
$AtomFeature{"HISND1"} = "Pos";
$AtomFeature{"HISNE2"} = "Pos";
$AtomFeature{"LYSNZ"} = "Pos";

$AtomFeature{"ASPOD1"} = "Neg";
$AtomFeature{"ASPOD2"} = "Neg";
$AtomFeature{"GLUOE1"} = "Neg";
$AtomFeature{"GLUOE2"} = "Neg";

$AtomFeature{"ALAO"} = "Acc";
$AtomFeature{"ARGO"} = "Acc";
$AtomFeature{"ASNO"} = "Acc";
$AtomFeature{"ASNOD1"} = "Acc";
$AtomFeature{"ASPO"} = "Acc";
$AtomFeature{"ASPOD1"} = "Acc";
$AtomFeature{"ASPOD2"} = "Acc";
$AtomFeature{"CYSO"} = "Acc";
$AtomFeature{"GLNO"} = "Acc";
$AtomFeature{"GLNOE1"} = "Acc";
$AtomFeature{"GLUO"} = "Acc";
$AtomFeature{"GLUOE1"} = "Acc";
$AtomFeature{"GLUOE2"} = "Acc";
$AtomFeature{"GLYO"} = "Acc";
$AtomFeature{"HISO"} = "Acc";
$AtomFeature{"ILEO"} = "Acc";
$AtomFeature{"LEUO"} = "Acc";
$AtomFeature{"LYSO"} = "Acc";
$AtomFeature{"METO"} = "Acc";
$AtomFeature{"PHEO"} = "Acc";
$AtomFeature{"PROO"} = "Acc";
$AtomFeature{"SERO"} = "Acc";
$AtomFeature{"THRO"} = "Acc";
$AtomFeature{"TRPO"} = "Acc";
$AtomFeature{"TYRO"} = "Acc";
$AtomFeature{"VALO"} = "Acc";

$AtomFeature{"ALAN"} = "Don";
$AtomFeature{"ARGN"} = "Don";
$AtomFeature{"ARGNE"} = "Don";
$AtomFeature{"ARGNH1"} = "Don";
$AtomFeature{"ARGNH2"} = "Don";
$AtomFeature{"ASNN"} = "Don";
$AtomFeature{"ASNND2"} = "Don";
$AtomFeature{"ASNOD1"} = "Don";
$AtomFeature{"ASPN"} = "Don";
$AtomFeature{"CYSN"} = "Don";
$AtomFeature{"GLNN"} = "Don";
$AtomFeature{"GLNNE2"} = "Don";
$AtomFeature{"GLUN"} = "Don";
$AtomFeature{"GLYN"} = "Don";
$AtomFeature{"HISN"} = "Don";
$AtomFeature{"HISND1"} = "Don";
$AtomFeature{"HISNE2"} = "Don";
$AtomFeature{"ILEN"} = "Don";
$AtomFeature{"LEUN"} = "Don";
$AtomFeature{"LYSN"} = "Don";
$AtomFeature{"LYSNZ"} = "Don";
$AtomFeature{"METN"} = "Don";
$AtomFeature{"PHEN"} = "Don";
$AtomFeature{"PRON"} = "Don";
$AtomFeature{"SERN"} = "Don";
$AtomFeature{"SEROG"} = "Don";
$AtomFeature{"THRN"} = "Don";
$AtomFeature{"THROG1"} = "Don";
$AtomFeature{"TRPN"} = "Don";
$AtomFeature{"TRPNE1"} = "Don";
$AtomFeature{"TYRN"} = "Don";
$AtomFeature{"TYROH"} = "Don";
$AtomFeature{"VALN"} = "Don"; 

$AtomFeature{"HISCG"} =  "Aro";
$AtomFeature{"HISND1"} =  "Aro";
$AtomFeature{"HISCD2"} =  "Aro";
$AtomFeature{"HISCE1"} =  "Aro";
$AtomFeature{"HISNE2"} =  "Aro";
$AtomFeature{"PHECG"} =  "Aro";
$AtomFeature{"PHECD1"} =  "Aro";
$AtomFeature{"PHECD2"} =  "Aro";
$AtomFeature{"PHECE1"} =  "Aro";
$AtomFeature{"PHECE2"} =  "Aro";
$AtomFeature{"PHECZ"} =  "Aro";
$AtomFeature{"TRPCG"} =  "Aro";
$AtomFeature{"TRPCD1"} =  "Aro";
$AtomFeature{"TRPCD2"} =  "Aro";
$AtomFeature{"TRPNE1"} =  "Aro";
$AtomFeature{"TRPCE2"} =  "Aro";
$AtomFeature{"TRPCE3"} =  "Aro";
$AtomFeature{"TRPCZ2"} =  "Aro";
$AtomFeature{"TRPCZ3"} =  "Aro";
$AtomFeature{"TRPCH2"} =  "Aro";
$AtomFeature{"TYRCD1"} =  "Aro";
$AtomFeature{"TYRCD2"} =  "Aro";
$AtomFeature{"TYRCE1"} =  "Aro";
$AtomFeature{"TYRCE2"} =  "Aro";
$AtomFeature{"TYRCG"} =  "Aro";
$AtomFeature{"TYRCZ"} =  "Aro";

$AtomFeature{"CYSS"} = "Sul";
$AtomFeature{"METSD"} = "Sul";


# ____________________________________________________________________________________________________________________
my $inFile = $ARGV[0];		# Input file
my $outFile = $ARGV[1];		# Output file
my $cutoff_step = $ARGV[2];	# Cutoff step
my $cutoff_limit = $ARGV[3];	# Cutoff limit

my $signType = $ARGV[4];	# Signature type:
							# 0:aCSM 
							# 1:aCSM-HP		{Hydro,Polar}
							# 2:aCSM-ALL		{Hydro,Pos,Neg,Acc,Don,Aro,Sul,Neutral}

my %acsm_all;
$acsm_all{"Hydro"} = 1;
$acsm_all{"Pos"} = 1;
$acsm_all{"Neg"} = 1;
$acsm_all{"Acc"} = 1;
$acsm_all{"Don"} = 1;
$acsm_all{"Aro"} = 1;
$acsm_all{"Sul"} = 1;
$acsm_all{"Neutral"} = 1;

if(scalar(@ARGV) < 5){
	print "\n###############################################################################################";
	print "\nSINTAX:\n\tperl aCSM.pl <infile> <outfile> <cutoff_step> ";
	print "<cutoff_limit> <signType{0,1,2}>\n\nEXAMPLE:\n\tperl aCSM.pl infile outfile 0.2 20.0 1\n\n";
	print "<infile> has the path to the pockets in PDB format.\n";
	print "<signType{0,1,2}> choose 0 for 'aCSM', 1 for 'aCSM-HP' or 2 for 'aCSM-ALL'.\n\n";	
	print "This example command will generate aCSM-HP signatures considering \na cutoff range of 0-20A, with a 0.2A step.\n";
	print "#################################################################################################";
	print "\n\n";
	exit;
}

if($cutoff_step < 0.001){
	print "\nError: Cutoff Step too small.\nChoose another please.\n\n";
	exit;
}

if($cutoff_limit > 1000){
	print "\nError: Cutoff Limit too large.\nChoose another please.\n\n";
	exit;
}

# ____________________________________________________________________________________________________________________
open(INFILE_1,"<$inFile") or die "$!Erro ao abrir: $inFile\n";
my @IDs = <INFILE_1>;
close INFILE_1;

open(OUTFILE,">$outFile") or die "$!Error: $outFile\n";

foreach my $pdbID (@IDs){
	chomp($pdbID);

	###################################
	print "Processing: $pdbID\n";
	###################################

	open(PDB,"<$pdbID") or die "$!Erro ao abrir: $pdbID\n";
	my @pdb = <PDB>;
	close PDB;

	my $max_index = 0;		my $k = 0;		my $temp = 0;
	my $sign = "";

	#--------------------------------------------
	# aCSM
	my $edgeCount = 0;
	# aCSM-HP
	my $edgeCountHH = 0;
	my $edgeCountPP = 0;
	my $edgeCountHP = 0;
	# aCSM-ALL
	my %edgeCountTipo2;

	my @keys = sort(keys(%acsm_all));
	foreach my $key1 (@keys){
		foreach my $key2 (@keys){
			$edgeCountTipo2{$key1.":".$key2} = 0;
			if($key1 ne $key2){
				$edgeCountTipo2{$key2.":".$key1} = 0;
			}
		}
	}
	#--------------------------------------------
	my $atom_selected = 0;
	my @coord_x;
	my @coord_y;
	my @coord_z;
	my @atom_index;
	my @res_name;
	my @atom_name;

	# Select atom lines
	foreach my $line (@pdb){
		if($line =~ m/^ATOM/){
			my $temp1 = trim(substr($line,30,8));
			my $temp2 = trim(substr($line,38,8));
			my $temp3 = trim(substr($line,46,8));
			$coord_x[$k] = $temp1;
			$coord_y[$k] = $temp2;
			$coord_z[$k] = $temp3;
			$atom_index[$k] = trim(substr($line,6,5));

			$res_name[$k] = substr($line,17,3);
			$atom_name[$k] = trim(substr($line,12,4));

			if($atom_index[$k] > $max_index){
				$max_index = $atom_index[$k];
			}
			$k++;
		}
	}

	$atom_selected = $max_index;
	my @dist;
	@keys = sort(keys(%acsm_all));

	# Cutoff limits
	my $cutoff_1 = 0.0;
	my $cutoff_2 = $cutoff_limit;
	# Calcula contatos e obtem numero de arestas
	for(my $i=0;$i<$k;$i++){
		for(my $j=0;$j<$k;$j++){
			# Somente parte superior da matriz calculada evitando calculos redundantes
			if($i < $j){
				$dist[$i][$j] = distance($coord_x[$i],$coord_y[$i],$coord_z[$i],
								$coord_x[$j],$coord_y[$j],$coord_z[$j]);
				if($dist[$i][$j] >= $cutoff_1 and $dist[$i][$j] <= $cutoff_2){
					$edgeCount++;

					my $key1 = $res_name[$i].$atom_name[$i];
					my $key2 = $res_name[$j].$atom_name[$j];

					if($signType == 1){
						if($Hydrophobic{$key1} and $Hydrophobic{$key2}){
							$edgeCountHH++;
						}
						elsif(!($Hydrophobic{$key1}) and !($Hydrophobic{$key2})){
							$edgeCountPP++;
						}
						else{
							$edgeCountHP++;
						}
					}
					elsif($signType == 2){
						if(!$AtomFeature{$key1}){
							$AtomFeature{$key1} = "Neutral";
						}
						if(!$AtomFeature{$key2}){
							$AtomFeature{$key2} = "Neutral";
						}

						$edgeCountTipo2{$AtomFeature{$key1}.":".$AtomFeature{$key2}}++;
						
						if($AtomFeature{$key1} ne $AtomFeature{$key2}){
							$edgeCountTipo2{$AtomFeature{$key2}.":".$AtomFeature{$key1}}++;
						}
					}
				}
			}
		}
	}

	if($signType == 0){	#aCSM
		$sign .= $edgeCount.",";
	}
	elsif($signType == 1){	#aCSM-HP
		$sign .= $edgeCountHH.",";
		$sign .= $edgeCountPP.",";
		$sign .= $edgeCountHP.",";
	}
	elsif($signType == 2){	#aCSM-ALL
		@keys = sort(keys(%acsm_all));
		for($a=0; $a<scalar(@keys); $a++){
			for($b=$a; $b<scalar(@keys); $b++){
				$sign .= $edgeCountTipo2{$keys[$a].":".$keys[$b]}.",";
			}
		}
	}

	undef %edgeCountTipo2;

	for(my $cutoff_temp=($cutoff_limit-$cutoff_step);$cutoff_temp>=0;$cutoff_temp-=$cutoff_step){
		$edgeCount = 0;
		$edgeCountHH = 0;
		$edgeCountPP = 0;
		$edgeCountHP = 0;
		undef %edgeCountTipo2;

		my @keys = sort(keys(%acsm_all));
		foreach my $key1 (@keys){
			foreach my $key2 (@keys){
				$edgeCountTipo2{$key1.":".$key2} = 0;
				if($key1 ne $key2){
					$edgeCountTipo2{$key2.":".$key1} = 0;
				}
			}
		}
	
		my $cutoff_1 = 0.0;
		my $cutoff_2 = ($cutoff_temp + 0.0);
		if($cutoff_2 < 0.0){
			last;
		}
		if($cutoff_2 <= ($cutoff_limit-$cutoff_step)*1.000001){
			for(my $i=0;$i<$k;$i++){
				for(my $j=0;$j<$k;$j++){
					if($i < $j){
						if($dist[$i][$j] >= $cutoff_1 and $dist[$i][$j] <= $cutoff_2){
							$edgeCount++;

							my $key1 = $res_name[$i].$atom_name[$i];
							my $key2 = $res_name[$j].$atom_name[$j];

							if($signType == 1){
								if($Hydrophobic{$key1} and $Hydrophobic{$key2}){
									$edgeCountHH++;
								}
								elsif(!($Hydrophobic{$key1}) and !($Hydrophobic{$key2})){
									$edgeCountPP++;
								}
								else{
									$edgeCountHP++;
								}
							}
							elsif($signType == 2){
								if(!$AtomFeature{$key1}){
									$AtomFeature{$key1} = "Neutral";
								}
								if(!$AtomFeature{$key2}){
									$AtomFeature{$key2} = "Neutral";
								}

								$edgeCountTipo2{$AtomFeature{$key1}.":".$AtomFeature{$key2}}++;
								#print "\n1. ".$AtomFeature{$key1}.": ".$key1."\n";
								#print "2. ".$AtomFeature{$key1}.": ".$key1."\n";
								#print "Dist: ".$dist[$i][$j]." - ".$cutoff_1." x ".$cutoff_2;


								if($AtomFeature{$key1} ne $AtomFeature{$key2}){
									$edgeCountTipo2{$AtomFeature{$key2}.":".$AtomFeature{$key1}}++;
								}
							}
						}
					}
				}
			}

			if($signType == 0){
				$sign .= $edgeCount.",";
			}
			elsif($signType == 1){
				$sign .= $edgeCountHH.",";
				$sign .= $edgeCountPP.",";
				$sign .= $edgeCountHP.",";
			}
			elsif($signType == 2){
				@keys = sort(keys(%acsm_all));
				for($a=0; $a<scalar(@keys); $a++){
					for($b=$a; $b<scalar(@keys); $b++){
						$sign .= $edgeCountTipo2{$keys[$a].":".$keys[$b]}.",";
					}
				}
			}
		}
	}
	# Next ID
	chop($sign);
	printf OUTFILE "$sign\n";
}
close OUTFILE;

exit;
# ____________________________________________________________________________________________________________________
sub trim{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub distance{
	my($x1,$y1,$z1,$x2,$y2,$z2) = @_;
	my $distance;
	
	$distance = sqrt(($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2);
	return $distance;
}
