# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 16:47:34 2021

@author: Jefer.Cano, Ricardo.Ortiz
"""

import pandas as pd
import numpy as np


#%%Cargue de los dataset base 

# Load verbatim-DwCA - Elementos priorizados + elementos que identifiques de GBIF
dwv = pd.read_table("D:\\HardDisck_Device\\ADATA\\Data\\SIB_Py\\Data\\Cumaral\\verbatim.txt", encoding = "utf8", sep= "\t", quoting=3, 
                    usecols = ['gbifID','institutionCode','collectionID', 'collectionCode','basisOfRecord','type','recordedBy','individualCount','eventID','eventDate','continent', 'country','stateProvince',
                               'county','municipality','minimumElevationInMeters','decimalLatitude','decimalLongitude','geodeticDatum','scientificName','scientificNameID','kingdom','genus','infraspecificEpithet','taxonRank','occurrenceStatus'])

# Cargar corte trimestral con datos para la validación
dwc_co = pd.read_table("D:\\HardDisck_Device\\ADATA\\Data\\SIB_Py\\Data\\Trimestre\\dwc_co_2022T4.txt", encoding = "utf8", sep= "\t", quoting=3, 
                    usecols = ['gbifID','publishingOrgKey', 'stateProvinceValidation','Departamento-ubicacionCoordenada','countyValidation','Municipio-ubicacionCoordenada','ZonaMaritima', 'flagGEO',
                               'flagTAXO','repatriated','publishingCountry','datasetKey','species','decimalLatitude','decimalLongitude','issue'])

# Join the datasets through GBIFID
dwm=dwc_co.merge(dwv, on='gbifID', suffixes=('_ocu', '_ver'))

#STOP... Remove the dwi and dwo variables from the explorer
del [dwc_co, dwv]

#Remove integrated datadatasets if you only want to run the diagnostic for national publishing data
#remove eBird from datasets
dwm=dwm[dwm.datasetKey != '4fa7b334-ce0d-4e88-aaae-2e0c138d049e']
#remove Xeno-Canto from dataset
dwm=dwm[dwm.datasetKey != 'b1047888-ae52-4179-9dd5-5448ea342a24']
#remove iNat from dataset
dwm=dwm[dwm.datasetKey != '50c9509d-22c7-4a22-a47d-8c48425ef4a7']
#remove repatriated from datasets
dwm=dwm[dwm.publishingCountry == 'CO']

#%% CHECK POINT
#Guardado de seguridad: datasets unidos únicamente con los datos de publicación nacional
dwm.to_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm.csv')
#Cargue de seguridad: cargue del dataset guardado en la linea anterior
dwm=pd.read_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm.csv', index_col=0) 

#%%
#Cargue de colecciones actualizadas para validacion y flags de los elementos priorizados DwC

#Cargar archivo de colecciones validadas actualizado. Extraer el datasetKey y la columna de validación
colecciones = pd.read_table("D:\HardDisck_Device\ADATA\Data\SIB_Py\Data\Colecciones\Colecciones_20220831.csv", encoding = "utf8", sep= ";", usecols = ['key','ColeccionVerificada'])
#Cargar archivo de reporte mensual actualizado. Extraer el datasetKey e IPT de la publicación.
datasetCO = pd.read_table("D:\HardDisck_Device\ADATA\Data\SIB_Py\Data\Colecciones\datasetco_20220930.txt", encoding = "utf8", sep= "\t", usecols = ['key','IPT','type','DOI_URL','title', 'NombreCorto_Org','year','typeOrg' ])    
#Union de las columnas para validacion de los registros biológicos del dataset dwm
dwm=pd.merge(left=dwm, right=colecciones, how= 'left', left_on='datasetKey', right_on='key')
dwm=pd.merge(left=dwm, right=datasetCO, how= 'left', left_on='datasetKey', right_on='key')
#Renombra las columnas agregadas y elimina la columna key duplicada
dwm.rename(columns={'type_x': 'type'}, inplace=True)
dwm.rename(columns={'type_y': 'CoreType'}, inplace=True)
dwm.rename(columns={'key_x': 'key'}, inplace=True)
del dwm['key_y']

#Verificación de las columnas nuevas
list(dwm)
#Elimina las variables no usadas
del[colecciones]

#%%CHECK POINT
#Guardado de seguridad: dwm unido con los dataset de las colecciones actualizadas
dwm.to_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm_coleccionesCumaral.csv')
#Cargue de seguridad: cargue del dataset guardado en la linea anterior
dwm=pd.read_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm_colecciones.csv', index_col=0)

#%%
#GENERAR VALIDACIONES Y FLAGS DE COLECCIONES Y DOCUMENTACION
#IMPORTANTE EN CASO DE: """AttributeError: Can only use .str accessor with string values!""", ejecutar las siguientes lineas antes con el fin de convertir estas columnas a tipo string, debido a que si no hay elementos tipo string no generara el flag 'flagGEO_textoenCoordenadas'.
dwm['decimalLatitude_ver']=dwm['decimalLatitude_ver'].astype(str)
dwm['decimalLongitude_ver']=dwm['decimalLongitude_ver'].astype(str)

#Crear columna de validación de texto en coordenadas geográficas decimales
#NOTA: se deben usar las columnas decimalLatitude_ver y decimalLongitude_ver correspondientes al verbatim para validar la presencia de caracteres de texto
dwm['flagGEO_textoenCoordenadas'] = np.where((dwm['decimalLatitude_ver'].str.contains(',| |N|S|E|W|O|°')) | (dwm['decimalLongitude_ver'].str.contains(',| |N|S|E|W|O|°')), 'Coordenada inválida: contiene valores de texto', '')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagGEO_textoenCoordenadas']).gbifID.nunique().reset_index()

#Crear columna de validación de vocabularios controlados en basisOfRecord
dwm['flagBasisOfRecord'] = np.where((dwm['basisOfRecord'].str.contains('HumanObservation|PreservedSpecimen|MachineObservation|FossilSpecimen|MaterialSample|LivingSpecimen')), '', 'Documentación inválida del elemento basisOfRecord')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagBasisOfRecord']).gbifID.nunique().reset_index()

#Crear columna de validación de vocabularios controlados en type
dwm['flagType'] = np.where((dwm['type'].str.contains('Objeto físico|Evento|Sonido|Imagen estática|Imagen en movimiento')), '', 'Documentación inválida del elemento type')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagType']).gbifID.nunique().reset_index()

#Crear columna de validación de vocabularios controlados en collectionID
dwm['flagCollectionID'] = np.where((dwm['ColeccionVerificada'].notnull()) & (dwm['collectionID'].str.contains('RNC:|Registro Nacional de Colecciones:|Registro Nacional de Colecciones Biológicas:')), '', 'Documentación inválida del elemento colletionID')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagCollectionID']).gbifID.nunique().reset_index()

#Crear columna de validación de vocabularios controlados en scientificNameID con IPT SiBM
dwm['flagscientificNameID'] = np.where((dwm['IPT']=='sibm') & (dwm['scientificNameID'].isnull()), 'Datos marinos: scientificNameID no documentado', '')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagscientificNameID']).gbifID.nunique().reset_index()

#Crear columna de validación de documentación de occurrenceStatus en datos con IPT SiBM
dwm['flagOccurrenceStatus_SiBM'] = np.where((dwm['IPT']=='sibm') & (dwm['occurrenceStatus'].isnull()), 'Datos marinos: occurrenceStatus no documentado', '')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagOccurrenceStatus_SiBM']).gbifID.nunique().reset_index()

#Crear columna de validación de vocabularios controlados en occurrenceStatus
dwm['flagOccurrenceStatus'] = np.where((~dwm['occurrenceStatus'].str.contains('present|absent', na = False)) & (dwm['occurrenceStatus'].notnull()), 'Documentación inválida en occurrenceStatus', '')
#Validación rápida. No de registros válidos e inválidos
dwm.groupby(['flagOccurrenceStatus']).gbifID.nunique().reset_index()

#Verificación de las columnas nuevas
list(dwm)

#%%CHECK POINT
#Guardado de seguridad: verificación de flags de colecciones y documentación
dwm.to_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm_flagSIB.csv')
#Cargue de seguridad: cargue del dataset guardado en la linea anterior
dwm=pd.read_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm_flagSIB.csv', index_col=0)

#%%
#Cargue de seguridad: cargue del dwm de los datos de colombia suministrado por ricardo, NOTA: registra dos issues que ya no existen en los generados por GBIF ()
dwm=pd.read_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Data\Colombia\dwm.csv', index_col=0)
list(dwm)
#%%
#GENERAR REPORTE FINAL DE FLAGS E ISSUES
#Filter by conditional type (Registros biologicos, Eventos de muestreo)
dwm = dwm[dwm['CoreType'].isin(['Registros biologicos','Eventos de muestreo'])] 

#Concatenar los flags del SIB y los issue de GBIF en la columna flags     
dwm['flags']=dwm['flagTAXO'].astype(str)+ ';'+dwm['flagGEO'].astype(str)+ ';'+dwm['flagGEO_textoenCoordenadas'].astype(str)+ ';'+dwm['flagBasisOfRecord'].astype(str)+ ';'+dwm['flagType'].astype(str)+ ';'+dwm['flagCollectionID'].astype(str)+ ';'+dwm['flagscientificNameID'].astype(str)+ ';'+dwm['flagOccurrenceStatus_SiBM'].astype(str)+ ';'+dwm['flagOccurrenceStatus'].astype(str)+ ';'+dwm['issue'].astype(str)
#Dejar solo columnas de interes para el repote final
dwm = dwm[['gbifID','datasetKey','flags']] 

#Explode the list of flags for every record in a issue or flag by row
dwm_explode = dwm.assign(flags=dwm['flags'].str.split(';')).explode('flags')

#Count the number of flags and issues grouped by datasetKey and flags 
flags_issues = dwm_explode.groupby(['datasetKey','flags']).gbifID.nunique().reset_index()
#Organize the table by datasetKey to make it more readable
final_dataset = flags_issues.pivot(index ='datasetKey', columns ='flags', values = 'gbifID').reset_index()

#Elimina las variables no usadas
del[flags_issues]

#%%
# Guaradado de seguridad: Exportar el dwm_explode al cual le asigno los valores de peso posteriormente.
dwm_explode.to_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\dwm_explode.csv', sep="\t", encoding = "utf8")
#Cargue de seguridad: cargue del dataset guardado en la linea anterior.
dwm_explode = pd.read_table("F:\ADATA\Data\SIB_Py\Check\dwm_explode.csv", index_col=0, encoding = "utf8", sep= "\t", quoting=3)

#%%
#Asignacion de valores de peso a los issues & flags para la priorizacion de los dataset.
#Cargue del archivo con los valores estandar para cada issue & flag priorizado, NOTA: el encoding queda en ANSI porque al corvetir el archivo en UTF-8 no me reconoce las tildes y por tanto no evaluara sobre los issues reportados.
valores_issues = pd.read_table("D:\\HardDisck_Device\\ADATA\\Data\\SIB_Py\\Data\\Colombia\\pesos_issues_final.txt", encoding = "ANSI", sep= "\t", quoting=3)
#Union del dwm_explode con los issues de la descarga y el valor_issues con los valores a asignar para cada issue de esta manera asigna los valores a los issues reportados, NOTA: No es necesario asignar valores a aquellos nan o vacios. 
pesos_asignados = pd.merge(left=dwm_explode, right=valores_issues, right_on='Issues', left_on='flags', how='inner')

#%%
#Crear los df con las sumas de los pesos para cada dataset.
#Agrupa los pesos por los dataset para poder hacer el merge con el final_dataset
flags_sum_dataset = pesos_asignados.groupby('datasetKey')['Total'].sum().reset_index()

#Suma los pesos de agrupando por el gbifID no por dataset, NOTA: no se agrega en el reporte final pero sirve para verificar que se hayan sumado correctamente.
flags_sum = pesos_asignados.groupby(['gbifID','datasetKey'])['Total'].sum().to_frame().reset_index()
#Elimina las variables no usadas
del[flags_sum]
#%%
#Merge de la suma de pesos con el final_dataset, y finalmente con el datasetCO para la informacion ligada a cada dataset
#Filter by conditional type (Registros biologicos, Eventos de muestreo)
datasetCO = datasetCO[datasetCO['type'].isin(['Registros biologicos','Eventos de muestreo'])] 

#Merge de la suma de pesos y el final dataset en base a el datasetKey
dwmF=pd.merge(left=datasetCO, right=final_dataset, how='left', left_on='key', right_on='datasetKey').merge(flags_sum_dataset, on='datasetKey')

#%%
#Export the final result to csv file.
dwmF.to_csv('D:\HardDisck_Device\ADATA\Data\SIB_Py\Check\ReporteFinal_PriorizacionColombia_2.csv', sep="\t", encoding = "utf8")


dwm.groupby([dwm['key']=!'9850b85a-793f-429d-a437-2d5b46282acf']).gbifID.nunique().reset_index()
