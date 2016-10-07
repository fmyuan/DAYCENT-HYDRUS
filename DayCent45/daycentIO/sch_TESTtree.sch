1800          Starting year
2005          Last year
site_tree.100   Site file name
0             Labeling type
-1            Labeling year
-1.00         Microcosm
-1            CO2 Systems
-1            pH shift
-1            Soil warming
0             N input scalar option
0             OMAD scalar option
0             Climate scalar option
2             Initial system
              Initial crop
MCON          Initial tree

Year Month Option
1             Block #   Spin-up (trees)
1898          Last year
1             Repeats # years
1800          Output starting year
1             Output month
0.083         Output interval
F             Weather choice
TEST_Daily.wth       
   1    1 TREE MCON 
   1    1 TFST
   1  365 TLST
-999 -999 X
2             Block #   Spin-up (fire)
1899          Last year
1             Repeats # years
1899          Output starting year
1             Output month
0.083         Output interval
C             Weather choice
   1    1 TREE MCON 
   1    1 TFST
   1  255 FIRE WILD 
   1  255 TREM WILD 
   1  365 TLST
-999 -999 X
3             Block #   pre-run
2005          Last year
1             Repeats # years
1900          Output starting year
1             Output month
0.083         Output interval
C             Weather choice
   1    1 TREE MCON 
   1    1 TFST
   1  365 TLST
-999 -999 X