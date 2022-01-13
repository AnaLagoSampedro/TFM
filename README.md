# TFM
Trabajo Fin de Máster Bioinformática avanzada UPO

Realizado por Ana María Lago Sampedro

Curso 2020-2021

Resumen del trabajo:

El avance en la tecnología aplicada a la Genómica ha supuesto un incremento en la cantidad de información obtenida que la llevó a ser considerada una disciplina más del Big Data y un escenario ideal para la aplicación de técnicas Machine Learning (ML). Debido a la naturaleza multigénica de enfermedades complejas como la Diabetes tipo 2 (DM2), es interesante enfocar su estudio como consecuencia de fallos en complejos módulos funcionales (o rutas de señalización metabólicas) causados por alteraciones en la actividad de los genes codificantes de proteínas interconectadas que los constituyen. Empleando novedosos modelos mecanicistas en combinación con estrategias ML se pueden explorar estos módulos a partir de datos de expresión génica para modelar el mapa mecanicista de la DM2 que señale aquellos circuitos (entendidos como sub-rutas) desregulados en la enfermedad y descubrir mediante aprendizaje automático, empleando grandes bases de datos de expresión génica de tejidos sanos, nuevos enfoques terapéuticos al aprender y predecir el impacto de genes diana de fármacos conocidos (KDTs) aprobados sobre esos circuitos y aplicar este conocimiento para el reposicionamiento de fármacos.



Aclaraciones:

En este repositorio se encuentran todos los scripts empleados para la realización de este trabajo, así como las tablas suplementarias y figuras indicadas en el trabajo escrito obtenidas. 

También se incluyen como anexo los scripts empleados para el análisis de expresión diferencial de los datos empleados para la realización de este trabajo y la tabla de resultados. Para el script Script_KDTs_relevantes.R se añaden las tablas empleadas con los resultados del ML; performance_stability_results_symbols.tsv, shap_summary_symbols.tsv y shap_selection_symbols.tsv

También se añaden la lista de rutas fisiológicas KEGG empleadas para la elaboración del mapa (physiological_paths.tsv) y los datos descargados de DrugBank (drugbank_drug-bindings_v5.1.8.tsv). 

Para obtener las matrices de datos empleadas obtenidas de los repositorio GEO y GTEx, en el propio trabajo se indica de donde descargarlas.
