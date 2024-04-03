# ROOT-analyser
Plusieurs version du programme pour l'analyse de données issues de LHcb en 2016,2017,2018,2019, DD ou LL et la combinaisons de plusieurs années. Pour plus de détails, le rapport est disponible (en anglais). Programmes par Bonafous Erwan & Couriol Florian


Quelques résultats: 
![Results2017KLL](https://github.com/DanmakuGecko/ROOT-analyser/assets/72706524/033ae54b-c38f-436c-89c6-d496667a799a)
![Screenshot from 2024-03-27 10-28-43](https://github.com/DanmakuGecko/ROOT-analyser/assets/72706524/75a42398-6d19-4448-bbfc-ef928e45f643)
![Screenshot from 2024-03-27 10-37-05](https://github.com/DanmakuGecko/ROOT-analyser/assets/72706524/f19b90b9-d697-4756-abb3-157e2833ebed)

cette fonction exploite une loi de puissance pour capturer les effets de perte d'énergie radiative. Cette forme s'aligne bien avec la distribution des énergies des particules lorsque la perte d'énergie radiative est un facteur prédominant. En combinant cette loi de puissance avec une distribution normale, nous obtenons une représentation plus précise des données expérimentales.

Normalement, pour décrire le bruit de fond, nous utilisons une exponentielle, mais dans nos données, nous avons observé que nous pouvons approximer l'exponentielle par une droite car nous avons peu de candidats de bruit de fond combinatoire. Et plus précisément par...

Pour modéliser correctement nos données, nous avons dû ajouter le fait que certains désintégrations de particules à quatre corps pourraient être détectées comme des désintégrations à trois corps. Les désintégrations peuvent être enregistrées comme une désintégration à trois corps car le π0π0 ne sera pas détecté, seuls les trois autres particules le seront. Pour intégrer cette désintégration dans notre modèle, nous avons ajouté deux "cristal balls" (pour et ), (il reste à dire que nous avons soustrait la masse).
