function [data, auxData, metaData, txtData, weights] = mydata_Ostrea_edulis

% http://www.debtheory.org/wiki/index.php?title=Mydata_file#Metadata

%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Ostreidae';
metaData.species    = 'Ostrea_edulis'; 
metaData.species_en = 'European Flat Oyster'; 
metaData.T_typical  = C2K(20); % K, body temp
metaData.data_0     = {'Lb', 'Ls'}; 
metaData.data_1     = {'t-L_F', 't-L_T', 't-JX_F', 't-JX_T', 't-Wd' };  

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011 http://www.debtheory.org/wiki/index.php?title=Completeness

metaData.author_1   = {'Brecht Stechele'};
metaData.date_subm_1 = [2020/02/10];
metaData.email_1    = {'Brecht.Stechele@ugent.be'};
metaData.address_1  = {'Ghent University Campus Coupure, Coupure links 653, 9000 Ghent, Belgium'};   

metaData.author_2  = {'Chenrui Yang'};
metaData.date_subm_2 = [2020/03/10]; 
metaData.email_2 = {'chenrui.yang@ugent.be'}; 
metaData.address_2 = {'Ghent University Campus Coupure, Coupure links 653, 9000 Ghent, Belgium'};
%% set data
% zero-variate data
data.Lb  = 0.014;      units.Lb  = 'cm';         label.Lb  = 'physical length at birth';           bibkey.Lb  = 'Helmer2019';  
data.Lr = 0.0183;      units.Lr  = 'cm';         label.Lr  = 'physical length at release';         bibkey.Lr  = 'Robert2017';  
data.Lv = 0.0300;   units.Lv  = 'cm';      label.Lv  = 'physical length at metamorphsis'; bibkey.Lv  = 'Robert2017';  
data.Ls = 0.0340;      units.Ls  = 'cm';         label.Ls  = 'physical length at settlement';      bibkey.Ls  = 'Naas1990';  
data.Li = 23;           units.Li = 'cm';    label.Li = 'maximul length';                 bibkey.Li  = 'Cole1941';

data.Wwp = 1.062;       units.Wwp = 'g';    label.Wwp = 'wet weight at puberty';         bibkey.Wwp = 'Elbe2013';  
data.Wdi = 35.8;        units.Wdi = 'g';    label.Wdi = 'ultimate dry weight';           bibkey.Wdi = 'Gals1964';
  comment.Wdi = 'based on the Crassostrea_virginica files';
%% uni-variate data

% Saurel 2003
% t-L
data.tL_Sau = [...
1       5.77
366     6.65
731     8.02
1096    7.84
1461    8.14
1826    8.77
2191    9.14
2556    9.38
2921    9.46
3286    10.29
3651    10.01
4016    9.63
4381    10.47
5111    11.09
5476    10.81
5841    11.26
6206    11.49
6571    11.5
7301    11.6
7666    12.19];
units.tL_Sau = {'d', 'cm'}; label.tL_Sau = {'time since start', 'total length'}; 
temp.tL_Sau    = C2K(12.4);  units.temp.tL_Sau = 'K';  label.temp.tL_Sau = 'temperature';
bibkey.tL_Sau = 'Saurel2003';

%% Robert 2017
%time vs ingestion

data.JX_30 = [...
0.00	28362.24
1.00	29903.66
2.01	42543.35
3.01	68131.02
4.00	101425.82
5.00	135645.47
5.99	52408.48
7.00	33294.80
8.00	27129.09
9.02	68131.02
10.0	72447.01
11.0	83853.56];
units.JX_30   = {'d', 'cell/d'};  label.JX_30 = {'time since release', 'number of ingested cell'};  
temp.JX_30   = C2K(30);  units.temp.JX_30 = 'K'; label.temp.JX_30 = 'temperature';
bibkey.JX_30 = 'Robert2017';

data.JX_25 = [...
0.00	20655.11
1.00	26204.24
2.01	44084.78
3.01	52100.19
4.00	71522.16
5.00	103892.10
6.01	73371.87
7.01	50867.05
8.00	36994.22
8.99	41001.93
9.99	45626.20
11.0	45934.49
12.0	95260.12    ];
units.JX_25   = {'d', 'cell/d'};  label.JX_25 = {'time since release', 'number of ingested cell'};  
temp.JX_25   = C2K(25);  units.temp.JX_25 = 'K'; label.temp.JX_25 = 'temperature';
bibkey.JX_25 = 'Robert2017';

data.JX_20 = [...
 0.01	9556.84
1.00	17263.97
1.99	25279.38
3.01	30211.95
4.00	36685.93
5.00	39768.79
5.99	58882.47
7.00	36685.93
8.00	39460.50
9.00	51483.62
10.0	43776.49
11.0	19113.68
12.0	13564.55
12.9	20346.82
14.0	10789.98
15.0	23737.96  ];
units.JX_20   = {'d', 'cell/d'};  label.JX_20 = {'time since release', 'number of ingested cell'};  
temp.JX_20   = C2K(20);  units.temp.JX_20 = 'K'; label.temp.JX_20 = 'temperature';
bibkey.JX_20 = 'Robert2017';

data.JX_15 = [...
0.01	5240.85
0.98	6782.27
2.01	5857.42
2.98	7398.84
4.00	9865.13
4.99	17880.54
6.01	21579.96
7.00	25279.38
8.00	20038.54
8.99	16955.68
10.0	15105.97
11.0	10173.41
12.0	10173.41
12.9	1849.71
14.0	616.57
15.0	1541.43 ] ;
units.JX_15   = {'d', 'cell/d'};  label.JX_15 = {'time since release', 'number of ingested cell'};  
temp.JX_15   = C2K(15); units.temp.JX_15 = 'K'; label.temp.JX_15 = 'temperature';
bibkey.JX_15 = 'Robert2017';
%% Roberta 2017
% t-JX_f
JX_f=[...
0.1     14225.79	21812.55	20614.31	15521.78
1   	20846.64	21795.46	23292.88	19997.72
1.97	19781.2     25671.89	26870.13	26271.39
2.96	19414.54	34042.11	36738.15	43328.47
3.94	25337.06	45707.49	55293.41	57390.34
4.96	20227.77	49584.3     54976.38	60368.47
5.98	22007.29	39681.34	42377.38	50165.95
6.96	21990.19	42659.47	46553.76	41161.67
7.98	18377.62	36650.8     47734.53	38448.16 
8.99	12368.95	22853.55	35135.53	32739.04
9.96	8457.57     27029.93	19540.92	36316.3
11.01	4245.5      49779.04	35999.27	57567.6
11.99	417         51259.37	50659.87	76721.99];
%500 cells/uL
data.JX_f500 = JX_f(:,[1,2]);
units.JX_f500 = {'d', 'cell/d'};  label.JX_f500 = {'time since release', 'length'};  
temp.JX_f500   = C2K(25);  units.temp.JX_f500 = 'K'; label.temp.JX_f500 = 'temperature';
bibkey.JX_f500 ='Robert2017';
%1500 cells/uL
data.JX_f1500 = JX_f(:,[1,3]);
units.JX_f1500 = {'d', 'cell/d'};  label.JX_f1500 = {'time since release', 'length'};  
temp.JX_f1500   = C2K(25);  units.temp.JX_f1500 = 'K'; label.temp.JX_f1500 = 'temperature';
bibkey.JX_f1500 ='Robert2017';
%2500 cells/uL
data.JX_f2500 = JX_f(:,[1,4]);
units.JX_f2500 = {'d', 'cell/d'};  label.JX_f2500 = {'time since release', 'length'};  
temp.JX_f2500   = C2K(25);  units.temp.JX_f2500 = 'K'; label.temp.JX_f2500 = 'temperature';
bibkey.JX_f2500 ='Robert2017';
%3500 cells/uL
data.JX_f3500 = JX_f(:,[1,5]);
units.JX_f3500 = {'d', 'cell/d'};  label.JX_f3500 = {'time since release', 'length'};  
temp.JX_f3500   = C2K(25);  units.temp.JX_f3500 = 'K'; label.temp.JX_f3500 = 'temperature';
bibkey.JX_f3500 ='Robert2017';

%% Robert 2017 
% time vs length

data.tL_15  = [...
0.01	0.0183
2.02	0.0191
4.03	0.0201
6.03	0.0208
8.01	0.0229
10.0	0.0239
11.0	0.0245
12.0	0.0245
13.0	0.0248
15.0	0.0248];
units.tL_15   = {'d', 'cm'};  label.tL_15 = {'time since release', 'total length'};  
temp.tL_15    = C2K(15); units.temp.tL_15 = 'K'; label.temp.tL_15 = 'temperature';
bibkey.tL_15 ='Robert2017';

data.tL_20 = [...
0.02    0.0183
2.02    0.0202
4.01	0.0230
6.03	0.0247
7.01	0.0262
8.03	0.0267
9.01	0.0279
10.01	0.0282];  
units.tL_20   = {'d', 'cm'};  label.tL_20 = {'time since release', 'total length'};  
temp.tL_20    = C2K(20);  units.temp.tL_20 = 'K'; label.temp.tL_20 = 'temperature';
bibkey.tL_20 = 'Robert2017';

data.tL_25 = [...
0.02	0.018358
2.00	0.020887
4.01	0.025075
6.03	0.029453
7.01	0.031302
8.03	0.032245 ];  % cm, total length at f and T
units.tL_25   = {'d', 'cm'};  label.tL_25 = {'time since release', 'total length'};  
temp.tL_25   = C2K(25);  units.temp.tL_25 = 'K'; label.temp.tL_25 = 'temperature';
bibkey.tL_25 = 'Robert2017';

data.tL_30 = [...
0.01	0.018321
2.02	0.021189
4.03	0.026547
5.03	0.029943
6.03	0.031717
8.01	0.032736 ];  % cm, total length at f and T
units.tL_30   = {'d', 'cm'};  label.tL_30 = {'time since release', 'total length'};  
temp.tL_30    = C2K(30);  units.temp.tL_30 = 'K'; label.temp.tL_30 = 'temperature';
bibkey.tL_30 = 'Robert2017';

%% Robert 2017 
% t-L_f
% food levels  = 500, 1500, 2500, 3500 cell/uL

tL_food =[...
0.1	0.017301	0.017301	0.017272	0.017301	0.017301
1	0.017893	0.018011	0.017863	0.017863	0.017331
2	0.019549	0.019372	0.019697	0.019519	0.017686
3	0.021057	0.021294	0.021294	0.020969	0.017686
4	0.021767	0.023364	0.023128	0.023482	0.017686
5	0.023216	0.025346	0.025375	0.025582	0.017686
6	0.023719	0.026706	0.026233	0.027091	0.017686
7	0.024961	0.028007	0.027948	0.028481	0.017804
8	0.025523	0.028333	0.028747	0.028865	0.017804
9	0.026203	0.028835	0.029043	0.02993     0.017804
10	0.026055	0.028185	0.028983	0.029013 	0.017774];

%500 cells/uL
data.tL_f500 = tL_food(:,[1,2]);
units.tL_f500   = {'d', 'cm'};  label.tL_f500 = {'time since release', 'length'};  
temp.tL_f500   = C2K(25);  units.temp.tL_f500 = 'K'; label.temp.tL_f500 = 'temperature';
bibkey.tL_f500 = 'Robert2017';

%1500 cells/uL
data.tL_f1500 = tL_food(:,[1,3]);
units.tL_f1500   = {'d', 'cm'};  label.tL_f1500 = {'time since release', 'length'};  
temp.tL_f1500   = C2K(25);  units.temp.tL_f1500 = 'K'; label.temp.tL_f1500 = 'temperature';
bibkey.tL_f1500 = 'Robert2017';

%2500 cells/uL
data.tL_f2500 = tL_food(:,[1,4]);
units.tL_f2500   = {'d', 'cm'};  label.tL_f2500 = {'time since release', 'length'};  
temp.tL_f2500   = C2K(25);  units.temp.tL_f2500 = 'K'; label.temp.tL_f2500 = 'temperature';
bibkey.tL_f2500 = 'Robert2017';

%3500 cells/uL
data.tL_f3500 = tL_food(:,[1,5]);
units.tL_f3500   = {'d', 'cm'};  label.tL_f3500 = {'time since release', 'length'};  
temp.tL_f3500   = C2K(25);  units.temp.tL_f3500 = 'K'; label.temp.tL_f3500 = 'temperature';
bibkey.tL_f3500 = 'Robert2017';

% %unfed level considered around 70 but food level is too low to make the model run at least the food level has to be 180
% data.tL_fl5 = Data.tL_fL(:,[1,6]);
% units.tL_fl5   = {'d', 'cm'};  
% label.tL_fl5 = {'time since release', 'length'};  
% temp.tL_fl5   = C2K(25);  
% units.temp.tL_fl5 = 'K'; 
% label.temp.tL_fl5 = 'temperature';
% bibkey.tL_fl5 = 'unfed level food t-L';

%% Larbarta 1999 
% t-Wd vs time, 
% LARVAE - to - POST-LARVAL GROWTH
% tL data at 16°C
data.tL_Laba = [0 7 11 17 18 20 23 27; 0.01736 0.02160 0.02490 0.0303 0.0345 0.0526 0.0727 0.0986 ]';
units.tL_Laba   = {'d', 'cm'};  label.tL_Laba = {'time since release', 'length'};  
temp.tL_Laba   = C2K(16);  units.temp.tL_Laba = 'K'; label.temp.tL_Laba = 'temperature';
bibkey.tL_Laba = 'Laba1999';

% tWd data at 16°C
data.tWd_Laba = [0 7 11 17 18 20 23 27; 1.7e-7 5.7e-7 1.15e-6 2.34e-6 2.26e-6 2.04e-6 4.65e-6 8.09e-6]'; 
units.tWd_Laba   = {'d', 'g'};  label.tWd_Laba = {'time since release', 'tissue dry weight'};  
temp.tWd_Laba   = C2K(16);  units.temp.tWd_Laba = 'K'; label.temp.tWd_Laba = 'temperature';
bibkey.tWd_Laba = 'Laba1999';

%% set weights for all real data
weights = setweights(data, []);
weights.Ls = 5* weights.Ls;
weights.Lr = 5* weights.Lr;
%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;

%% Group plots
set1 = {'JX_30','JX_25','JX_20', 'JX_15'}; comment1 = {'ingested cells at diff temperature from Robert 2017'};
set2 = {'JX_f500', 'JX_f1500', 'JX_f2500','JX_f3500'}; comment2 = {'ingested cells at diff food level from robert 2017'};
set3 = {'tL_15', 'tL_20', 'tL_25','tL_30'}; comment3 = {'growth at diff temp from robert 2017'};
set4 = {'tL_f500', 'tL_f1500','tL_f2500','tL_f3500'}; comment4 = {'growth at diff food level from Robert 2017'};
metaData.grp.sets = {set1, set2, set3,set4};
metaData.grp.comment = {comment1, comment2, comment3,comment4};

%% Facts
%F1 = '';
%metaData.bibkey.F1 = ''; 
%metaData.facts = struct('F1',F1);

%% Discussion points
D1 = '';
D2 = '';     
metaData.discussion = struct('D1', D1, 'D2', D2);

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/my_pet}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parametaers and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Helmer2019'; type = 'article'; bib = [ ... 
'author = {Helmer, Luke and Farrell, Paul and Hendy, Ian and Harding, Simon and Robertson, Morven and Preston, Joanne}, ' ...
'year = {2019}, ' ...
'title  = {Active management is required to turn the tide for depleted Ostrea edulis stocks from the effects of overfishing, disease and invasive species}, ' ...
'journal = {Peer}, ' ...
'volume = {7}, ' ...
'pages = {e6431}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Robert2017'; type = 'article'; bib = [ ... 
'author = {Robert, Rene and Vignier, Julien and Petton, Bruno}, ' ...
'year = {2017}, ' ...
'title  = {Influence of feeding regime and temperature on development and settlement of oyster Ostrea edulis (Linnaeus, 1758) larvae}, ' ...
'journal = {Aquaculture Research}, ' ...
'volume = {48}, ' ...
'pages = {4756-4773}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Labarta1999'; type = 'article'; bib = [ ... 
'author = {Labarta, Uxo and Fernandez-Reiriz, Maria Jose and Perez-Camacho, A}, ' ...
'year = {2017}, ' ...
'title  = {Energy, biochemical substrates and growth in the larval development, metamorphosis and postlarvae of Ostrea edulis (L.)}, ' ...
'journal = {Journal of Experimental Marine Biology and Ecology}, ' ...
'volume = {238}, ' ...
'pages = {225-242}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

bibkey = 'Naas1990'; type = 'article'; bib = [ ... 
'author = {Naas, Kjell E}, ' ...
'year = {2017}, ' ...
'title  = {A semi-intensive metahod for spat production of the European flat oyster (Ostrea edulis L.)}, ' ...
'journal = {Aquacultural engineering}, ' ...
'volume = {9}, ' ...
'pages = {447-451}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

bibkey = 'Saurel2003'; type = 'article'; bib = [ ... 
'author = {Saurel and Richardson}. ' ...
'year = {2003}. ' ...
'title  = {Age and growth analysis of native oyster bed (Ostrea edulis) in Wales}. ' ...
'journal = {CCW Contract Science Report No 549}. ' ...
'volume = {549}. ' ...
'pages = {-}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

%
