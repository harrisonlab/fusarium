.fastq to old style fastq

# test set
cat 125_250k_F_trim.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >125_250k_old_F_trim.fastq
cat 125_250k_R_trim.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >125_250k_old_R_trim.fastq

/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/newAssembly -force assembly/newbler/125_newbler
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler/ tmp/125_250k_old_F_trim.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler/ tmp/125_250k_old_Rtrim.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/runProject assembly/newbler/125_newbler

#full assembly (untrimmed)

cat tmp/125_S3_L001_R1_001.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >tmp/125_old_F.fastq
cat tmp/125_S3_L001_R2_001.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >tmp/125_old_R.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/newAssembly -force assembly/newbler/125_newbler_complete_large
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler_complete_large/ tmp/125_old_F.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler_complete_large/ tmp/125_old_R.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/runProject -cpu 16 -large assembly/newbler/125_newbler_complete_large
# readStatus
#         {
#                 numAlignedReads    = 11006659, 92.71%, 46.32%;
#                 numAlignedBases    = 3051857465, 92.28%, 86.02%;
#                 inferredReadError = 0.78%, 23859603;
# 
#                 numberAssembled = 10663828, 89.83%, 44.88%;
#                 numberPartial   = 342347, 2.88%, 1.44%;
#                 numberSingleton = 88960, 0.75%, 0.37%;
#                 numberRepeat    = 756296, 6.37%, 3.18%;
#                 numberOutlier   = 20127, 0.17%, 0.08%;
#                 numberTooShort  = 0, 0.00%, 0.00%;
#         }
# 
#         pairedReadStatus
#         {
#                 numberWithBothMapped   = 10709150;
#                 numberWithOneUnmapped  = 176500;
#                 numberMultiplyMapped   = 957422;
#                 numberWithBothUnmapped = 19288;
#			}
# scaffoldMetrics
#         {
#                 numberOfScaffolds   = 954;
#                 numberOfBases       = 50332019;
# 
#                 avgScaffoldSize     = 52758;
#                 N50ScaffoldSize     = 176120, 64;
#                 largestScaffoldSize = 1797911;
# 
#                 numberOfScaffoldContigs     = 1521;
#                 numberOfScaffoldContigBases = 50029421;
# 
#                 avgScaffoldContigSize       = 32892;
#                 N50ScaffoldContigSize       = 113967, 127;
#                 largestScaffoldContigSize   = 838859;
#		}
#         largeContigMetrics
#         {
#                 numberOfContigs   = 1824;
#                 numberOfBases     = 50352439;
# 
#                 avgContigSize     = 27605;
#                 N50ContigSize     = 113426;
#                 largestContigSize = 838859;
# 
#                 Q40PlusBases      = 50340986, 99.98%;
#                 Q39MinusBases     = 11453, 0.02%;
# 
#                 largeContigEndMetrics
#                 {
#                         NoEdges   = 551, 15.1%;
#                         OneEdge   = 2202, 60.4%;
#                         TwoEdges  = 498, 13.7%;
#                         ManyEdges = 78, 2.1%;
#                         LargeRep  = 319, 8.7%;
#                 }
#         }
# 
#         allContigMetrics
#         {
#                 numberOfContigs = 2899;
#                 numberOfBases   = 50625365;
#         }





#full assembly (trimmed)

cat tmp/125_S3_L001_R1_001.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >tmp/125_old_F.fastq
cat tmp/125_S3_L001_R2_001.fastq | awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s#0/%s (%s)\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7], substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' >tmp/125_old_R.fastq
fastq-mcf /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa tmp/125_old_F.fastq tmp/125_old_R.fastq -C 10000000 -u -k 20 -t 0.01 -q 30 -p 5
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/newAssembly -force assembly/newbler/125_newbler_complete_large_trim
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler_complete_large_trim/ tmp/125_old_F_trim.fastq
/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/addRun -lib 250_PE -p assembly/newbler/125_newbler_complete_large_trim/ tmp/125_old_R_trim.fastq

# readStatus
#         {
#                 numAlignedReads    = 10947885, 92.23%, 46.09%;
#                 numAlignedBases    = 2988056640, 91.79%, 89.12%;
#                 inferredReadError = 0.75%, 22381374;
# 
#                 numberAssembled = 10629633, 89.55%, 44.75%;
#                 numberPartial   = 317924, 2.68%, 1.34%;
#                 numberSingleton = 87185, 0.73%, 0.37%;
#                 numberRepeat    = 816218, 6.88%, 3.44%;
#                 numberOutlier   = 19719, 0.17%, 0.08%;
#                 numberTooShort  = 0, 0.00%, 0.00%;
#         }
# 
#         pairedReadStatus
#         {
#                 numberWithBothMapped   = 10671246;
#                 numberWithOneUnmapped  = 174856;
#                 numberMultiplyMapped   = 997058;
#                 numberWithBothUnmapped = 18044;
# 
#                 library
#                 {
#                         libraryName       = "250_PE";
#                         libraryNumPairs   = 5930602;
#                         numInSameScaffold = 3980265, 67.1%;
# 
#                         pairDistanceRangeUsed      = 337..1137;
#                         computedPairDistanceAvg    = 674.2;
#                         computedPairDistanceDev    = 231.6;
#                 library
#                 {
#                         libraryName       = "250_PE";
#                         libraryNumPairs   = 5931180;
#                         numInSameScaffold = 4069416, 68.6%;
# 
#                         pairDistanceRangeUsed      = 331..1115;
#                         computedPairDistanceAvg    = 663.2;
#                         computedPairDistanceDev    = 226.3;
#                 }
#         }
# scaffoldMetrics
#         {
#                 numberOfScaffolds   = 959;
#                 numberOfBases       = 50304984;
# 
#                 avgScaffoldSize     = 52455;
#                 N50ScaffoldSize     = 177723, 63;
#                 largestScaffoldSize = 1799507;
# 
#                 numberOfScaffoldContigs     = 1559;
#                 numberOfScaffoldContigBases = 49975322;
# 
#                 avgScaffoldContigSize       = 32056;
#                 N50ScaffoldContigSize       = 113324, 127;
#                 largestScaffoldContigSize   = 838865;
# 			}
# largeContigMetrics
#         {
#                 numberOfContigs   = 1897;
#                 numberOfBases     = 50319660;
# 
#                 avgContigSize     = 26525;
#                 N50ContigSize     = 110898;
#                 largestContigSize = 838865;
# 
#                 Q40PlusBases      = 50308231, 99.98%;
#                 Q39MinusBases     = 11429, 0.02%;
# 
#                 largeContigEndMetrics
#                 {
#                         NoEdges   = 576, 15.2%;
#                         OneEdge   = 2249, 59.3%;
#                         TwoEdges  = 527, 13.9%;
#                         ManyEdges = 83, 2.2%;
#                         LargeRep  = 359, 9.5%;
#                 }
#         }
# 
#         allContigMetrics
#         {
#                 numberOfContigs = 3007;
#                 numberOfBases   = 50593093;
#         }






/home/armita/prog/DataAnalysis_2.9_All/gsnewbler_2.9-2_amd64/opt/454/apps/mapper/bin/runProject -cpu 16 -large assembly/newbler/125_newbler_complete_large_trim