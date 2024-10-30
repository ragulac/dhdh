#!/usr/bin/python3
# -*- coding: latin-1 -*-

import glob,os
from datetime import datetime,timedelta
from math import sqrt
import numpy as np
from netCDF4 import Dataset, MFDataset
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
dutc=5

def user():
#        ruta='mmif2022-03'
#        ruta='datos2022'
#        ruta='septiembre'
#        ruta='datos_7_diciembre'
        ifecha=datetime(2023,1,1,0)
        efecha=datetime(2023,12,31,23)	
        dtlocal=-5

        namelist={}
#        namelist['ruta']=ruta
        namelist['ifecha']=ifecha
        namelist['efecha']=efecha
        namelist['dtlocal']=dtlocal
#        namelist['lista']='prueb2_abril.inp'
#        namelist['lista']='datos07ddfi.inp'
        namelist['lista']='mmif_2023.inp'
#        namelist['Estad']='logs/PrecDia_abril2022_prueba2.csv'
        namelist['Estad']='logs/PrecDia_2023.csv'


        return namelist


#---------------Funciones de cálculo---------------
def humrel(t2,q2,psfc):

    """
       input:
       Campos temperatura t(Celsius), razon de mezcla qv(kg/kg), presion(Pa)
       (np.arrays) 1D,2D,...

       return:
       np.array 1D,2D,... de la humedad relativa(%)
    """
    e = (psfc * q2) / (0.622 + q2)
    esat = 611.2 * np.exp((17.67 * (t2 - 273.15)) / (t2 - 29.65))
    return e / esat * 100.0
    
def err_abs(a,b):
    errabs=np.zeros(len(a))
    for i in range(len(a)):
        if(b[i]!=0):
            errabs[i]=(a[i]-b[i])/b[i]
    return np.average(errabs)

def RMSE(a,b):
    from math import sqrt
    rmse=0
    for i in range(len(a)):
    	rmse=rmse+(a[i]-b[i])*(a[i]-b[i])
    rmse=sqrt(rmse)/len(a)
    return rmse

def BIAS(a,b):
    bias=0
    for i in range(len(a)):
    	bias=bias+(a[i]-b[i])
    bias=bias/len(a)
    return bias

def angle_wind(U10,V10):
#calcular el angulo en grados del vector de componentes U10,V10
    pi=3.141592653589793
    radeg=180./pi
    alpha = np.arctan2(U10, V10)*radeg
    return alpha

def norm_wind(U10,V10):
#calcular la norma del vector con componentes U10, V10 
    norm=np.sqrt(U10**2+V10**2)
    return norm

def ll_xy(lat,lon,lat0,lon0):
# Punto de grilla mas cercano a (lat,lon)
    import os
    from scipy import io
    lat0=float(lat0)
    lon0=float(lon0)
    #il <-> lat s-n = y
    for il in range(1,len(lat)):
        #jl <-> lon w-e = x
        for jl in range(1,len(lon[il])):
#            print(lat[il][jl],lon[il][jl],lat0,lon0)
            if lat[il][jl]>lat0 and lat[il-1][jl]<=lat0:
                if lon[il][jl]>lon0 and lon[il][jl-1]<=lon0:
                    dlat=lat[il][jl]-lat[il-1][jl]
                    dlon=lon[il][jl]-lon[il][jl-1]
                    nxmin=jl-1
                    nymin=il-1
                    if lat[il][jl]<lat0+dlat/2.: nymin=il
                    if lon[il][jl]<lon0+dlon/2.: nxmin=jl
#                    print( ('Punto WRF | Obs:'),nxmin,nymin,lon[nymin][nxmin],lat[nymin][nxmin],' | ',lon0,lat0)
                    break
    return nxmin,nymin,lat[nymin][nxmin],lon[nymin][nxmin]		

def ll_xy_np(lat,lon,lat0,lon0):
    #axis0=i=lat
    dlats2d=np.abs(lat-lat0)
    latp=np.min(dlats2d,axis=0) #x-array de lats mas cercanas para todas las lons
    ijmin=np.where(dlats2d==latp) 
    indices2d=list(zip(ijmin[0],ijmin[1])) #(i,j) de 1D lats mas cercanas 
    lons=np.array([lon[ij[0],ij[1]] for ij in indices2d])  # lons de (i,j)'s  
    dlons1d=np.abs(lons-lon0) 
    lonp=np.where(dlons1d==np.min(dlons1d)) # punto j
    ij=indices2d[lonp[0][0]] # coos (i,j) del punto
    #i=sw,j=we
    nymin=ij[0]
    nxmin=ij[1]
#    print( ('Punto WRF | Obs:'),nxmin,nymin,lon[nymin][nxmin],lat[nymin][nxmin],' | ',lon0,lat0)
    return nxmin,nymin,lat[nymin][nxmin],lon[nymin][nxmin]	


#---------------Funciones de lectura---------------

def aptos():
#Lee el archivo de aeropuertos y extrae codigo ICAO, nombre de aeropuerto, longitud y latitud

        listaptos='aeropuertos_24.csv'
#        listaptos='aeropuertos_precdia3.csv'
#        listaptos='aeropuertos_SYNOP.csv'
#        listaptos='aeropuertos_humrel.csv'
#        listaptos='Rubiales.csv'
#        listaptos='aeropuertos.csv'
#        listaptos='aeropuertos_3dvar.csv'
#        listaptos='aeropuertos_aj.bk'
#        listaptos='aeropuertos_27.csv'
#        listaptos='aeropuertos_aj.csv'
#        listaptos='ingeominas.csv'
        a=open(listaptos,'r')
        aptos=a.readlines()
        coos={}
        for linea in aptos:
            linea=linea.rstrip('\n')
            linea=linea.split(';')
            nome=linea[0]
            coos[nome]=[] #Codigo ICAO
            coos[nome].append(linea[2]) #Longitud
            coos[nome].append(linea[3]) #Latitud
            coos[nome].append(linea[1]) #Nombre aeropuerto
	
	#pprint(coos)

        return coos		

def DataSynop(sitios):

    DataAer={} 
    for apto in sitios.keys():
#      try:
        print(apto)
        DataObs={} 
        loc='/home/wrf/inidata/SYNOP/'+apto+'2023.csv'
        a=open(loc,'r')
#        print(loc)
        data=a.readlines()
        #dict ordenado por fechas
        odata={}
        fechas=[]
#        print(data)
        for line in data[1:]:
            linea=line.rstrip('/n').split(';')
            ff=datetime.strptime(linea[0],"%Y-%m-%d %H:00") #fechas en UTC
            odata[ff]=linea
            fechas.append(ff)
        fechas.sort()
        #inicializa dia prec
        iff=fechas[0]
        idiap=datetime(iff.year,iff.month,iff.day,(7+dutc)) #dia lluvia en UTC
        ediap=idiap+timedelta(hours=24)
        prec24=0
#       2021-03-27 12:00;600;7;40;4;12.5;11.4;93;752.9;9999;12;24
        lasttpr=idiap
        lastpr=0
        for ff in fechas:
            linea=odata[ff]
            if ff>ediap:  
                hoy=datetime(ediap.year,ediap.month,ediap.day)
                DataObs[hoy]=prec24
                idiap=datetime(ff.year,ff.month,ff.day,(7+dutc))
                ediap=idiap+timedelta(hours=24)
#                print('Diaprec',apto,hoy,ff,prec24)
                prec24=0
                lasttpr=idiap
                lastpr=0
            if '9999' not in linea[10]:
                prec=float(linea[10])
                tpr=float(linea[11])
                itpr=ff-timedelta(hours=tpr)
                if itpr >= idiap:
                    if itpr>=lasttpr:
                        prec24=prec24+prec
                        lasttpr=ff
                        lastpr=prec
                    else:
                        prec24=prec24+prec-lastpr
                        lastpr=prec
        DataAer[apto]=DataObs
    return DataAer


def DataOBS(sitios):
#Lee el archivo de datos observacionales de cada uno de los 'sitios' provenientes de ./METAR transformados mediante metar.py
#Modificar si se desea extraer otra variable o si cambia el formato de los archivos de los datos de observacion

        DataAer={} #Diccionario doble: el key es igual al de 'sitios' y su valor son los datos horarios de observacion
        for apto in sitios.keys():
#		loc='/home/nikolas/Ajust/METAR/11meses/'+apto+'_2016.csv'
#                '/home/nikolas/Ajust/METAR/11meses/Ingeominas-2015.csv'
                loc='/home/wrf/inidata/METAR/'+apto+'_2023.csv'
#                loc='/tmp/Prec_'+apto+'.csv'
                a=open(loc,'r')
                data=a.readlines()
                Data_OBS={}
                for x in range(1,len(data)):
                        linea=data[x].rstrip('\n')
                        linea=linea.split(',')
                        fecha=datetime.strptime(linea[5],'%Y-%m-%d %H:00') #fecha para el key del diccionario
                        if '9999' not in linea[6]:
                                Temperatura=float(linea[6]) #Extraer temperatura
                                Data_OBS[fecha]=Temperatura
#			else:
#				print '### CasoVariable no existente ###'
#				Data_OBS[fecha]=9999
                DataAer[apto]=Data_OBS
#	pprint(DataAer)
        return DataAer

def extractfromwrf(nx,ny,archivo,HoraTemp):

#    archivo=Dataset(nombre_archivo,'r')
    wrf_lat=np.copy(archivo['XLAT'][0,ny,nx].data)
    wrf_long=np.copy(archivo['XLONG'][0,ny,nx].data)
    ValTemp1=np.copy(archivo['U10'][0,ny,nx].data)
    ValTemp2=np.copy(archivo['V10'][0,ny,nx].data)
    ValRes=norm_wind(ValTemp1,ValTemp2)
    angulo=angle_wind(ValTemp1,ValTemp2)
    del ValTemp1
    del ValTemp2c
#    archivo.close()
 
    return ValRes, wrf_lat, wrf_long		



#---------------Funciones de comparacion de datos observados y datos WRF---------------

def eval(Dicc_WRF,Dicc_Obs,log):
#Esta funcion correlaciona los datos de los diccionarios de WRF y de Observacion para el valor 'apto'
#Retorna el par de arreglos utilizados en el calculo de estadisticos

	#Inicializar el Log (opcional)
#        nlog='Viento.2019_'+apto+'.csv'
#        log=open(nlog,'w')
        #Inicializar arreglos de retorno
        wrf_st=[]
        obs_st=[]
        fechas=[]
	#Ciclo sobre elementos del diccionario de WRF
        fwrf=[*Dicc_WRF]
        fwrf.sort()
        for wrfkey in fwrf: #for sobre los WRF.keys()
                fecha=wrfkey
		#Si la fecha se encuentra tambien en las observaciones
                if fecha in [*Dicc_Obs]:
		    #Agregar a los arreglos
                    wrf_st.append(Dicc_WRF[fecha]) #sin el [0] el objeto es una lista!
                    obs_st.append(Dicc_Obs[fecha]) 
                    fechas.append(fecha)
                    log.writelines(str(fecha)+';'+str(Dicc_WRF[fecha])+';'+str(Dicc_Obs[fecha])+'\n')

        return wrf_st,obs_st,fechas

def aeropuerto(coos,latw,lonw,var,apto):

    #Coordenadas geoespaciales a punto de grilla
    lat=float(coos[apto][0])
    lon=float(coos[apto][1])
    nx,ny,wlat,wlon=ll_xy_np(latw,lonw,lat,lon)
    dato=var[ny,nx]

    return (dato,apto,wlat,wlon)		

def main():
        import multiprocessing as mp

        namelist=user()
        lista=namelist['lista']
#        ruta=namelist['ruta']
        ifecha=namelist['ifecha']
        efecha=namelist['efecha']
                
        DWRF={}
        log={}
        wrf_total=[]
        obs_total=[]
        nn=mp.cpu_count()

	#Abrir tabla de informacion compilada
        Estadisticos=namelist['Estad']
        r=open(Estadisticos,'w')
        r.writelines('aeropuerto;latitud;longitud;#datos;latitud_WRF;longitud_WRF;PromWRF;PromObs;Pearson;Spearman;RMSE;BIAS\n')

	#nombre y coos[apto] 
        coos=aptos()
        for apto in [*coos]:
            DWRF[apto]={}
            nlog='logs/Prec.2023_'+apto+'.csv'
            log[apto]=open(nlog,'w')

	#Series Obs series[apto]
#        series=DataOBS(coos)
#       series=DataOBS[OACI][datetime(dia)]
        series=DataSynop(coos) # fechas en UTC

        #Series WRF
        wrfin=open(lista,'r')
        wrfnc=wrfin.readlines()
        wlat={}
        wlon={}
        for li in range(0,len(wrfnc)):
            nombre=wrfnc[li].rstrip('\n')
#            print(nombre)
            archivo=Dataset(nombre,'r')
            latw=archivo['XLAT'][0]
            lonw=archivo['XLONG'][0]
            rainc=archivo['RAINC'][0]
            rainnc=archivo['RAINNC'][0]
#            var=rainc+rainnc
            var=rainnc
            DiaTemp=archivo['Times'][0]
            HoraTemp=[str(i) for i in DiaTemp]
            HoraTemp=''.join(HoraTemp)
            HoraTemp=HoraTemp.replace("'b'",'').lstrip("b'").rstrip("'")
            #todo en hora local
            HoraTemp=datetime.strptime(HoraTemp,'%Y-%m-%d_%H:00:%S')
            HoraTemp=datetime(HoraTemp.year,HoraTemp.month,HoraTemp.day,HoraTemp.hour)
#            -timedelta(hours=dutc) #series WRF hora local
            print('FECHA: ',HoraTemp)
            if HoraTemp <= efecha:
                with mp.Pool(processes=nn) as pr:  
                    varout=pr.starmap(aeropuerto,[(coos,latw,lonw,var,apto) for apto in [*series]])
                for da in varout: 
                    DWRF[da[1]][HoraTemp]=da[0]                
                    wlat[da[1]]=da[2]
                    wlon[da[1]]=da[3]
            else:
                break
            archivo.close()

        for apto in [*coos]:
            Sobs=series[apto] #series en UTC 
            Swrf=DWRF[apto]   #serie en UTC
            Pwrf={}
            Pobs={}
            #precip diaria WRF:
            fechas=list(Sobs)
            fechas.sort()
            day=0
            for fecha in fechas[1:]:
              if fecha > ifecha:
                if fecha > efecha-timedelta(hours=24): break
                year=fecha.year
                mes=fecha.month
                dia=fecha.day
                hoy=datetime(year,mes,dia,0)
                #dia pluv = 7a.m. a 7a.m. <=> 12Zhoy-00Zhoy+23Zayer-12Zayer = 7hoy-19ayer+18ayer-7ayer, corrida: 00Z-23Z
                #                     7a.m                -      7p.m.ayer                +         6p.m.ayer                   -         7a.m.ayer
                #                   00Z+12                       00Z                                 23Zayer                    -         12Zayer
#                Pwrf[hoy]=Swrf[hoy+timedelta(hours=7)]-Swrf[hoy-timedelta(hours=dutc)]+Swrf[hoy-timedelta(hours=dutc+1)]-Swrf[hoy-timedelta(hours=12+dutc)]
                Pwrf[hoy]=Swrf[hoy+timedelta(hours=12)]-Swrf[hoy]+Swrf[hoy-timedelta(hours=1)]-Swrf[hoy-timedelta(hours=12)]
                #Por si hay obs horarias
                Pobs[hoy]=0.
                if day==dia:
                    Pobs[hoy]=Pobs[hoy]+Sobs[fecha]
                else:
                    Pobs[hoy]=Sobs[hoy]
                    day=dia
            wrf_st,obs_st,fecha=eval(Pwrf,Pobs,log[apto])
            wrf_total=wrf_total+wrf_st
            obs_total=obs_total+obs_st
            print(apto,len(wrf_st),len(obs_st))
            if len(wrf_st)>2:
                lat=float(coos[apto][1])
                lon=float(coos[apto][0])
                correl=pearsonr(obs_st,wrf_st)[0]
                spear=spearmanr(wrf_st,obs_st)[0]
                r.writelines(apto+';'+str(lat)+';'+str(lon)+';'+str(len(wrf_st))+';'+str(wlat[apto])+';'+str(wlon[apto])+';'+str(sum(wrf_st)/len(wrf_st))+';'+str(sum(obs_st)/len(obs_st))+';'+str(correl)+';'+str(spear)+';'+str(RMSE(wrf_st,obs_st))+';'+str(BIAS(wrf_st,obs_st))+'\n')
#                plt.figure()
#                plt.plot(wrf_st, label='wrf')
#                plt.plot(obs_st, label='obs')
#                plt.legend()
#                plt.savefig(apto+'.png')
#                plt.close()
                print('wrf_st',len(wrf_st))
                print('obs_st',len(obs_st))
                print( 'Estacion',apto)
                print( 'No. datos',len(wrf_st),len(obs_st),sum(wrf_st),sum(obs_st))
                print( 'promedio',sum(wrf_st)/float(len(wrf_st)),sum(obs_st)/float(len(obs_st)))
                print( 'Pearson', pearsonr(obs_st,wrf_st))
                print( 'Spearman',spearmanr(wrf_st,obs_st))
                print( 'RMSE', RMSE(wrf_st,obs_st) )
                print( 'BIAS', BIAS(wrf_st,obs_st),'\n')
#                ST(wrf_st,obs_st,fecha,apto)
#            else:
#                print( '### Aeropuerto Fallido ### sin datos correlacionados','\n')

        r.writelines('\n')
        r.writelines( 'PEARSON TOTAL:'+str(pearsonr(obs_total,wrf_total)[0]))
        r.writelines('\n')
        r.writelines( 'SPEARMAN TOTAL:'+str(spearmanr(wrf_total,obs_total)[0]))
        r.writelines('\n')
        r.writelines( 'RMSE TOTAL'+str(RMSE(wrf_total,obs_total)))
        r.writelines('\n')
        r.writelines( 'BIAS TOTAL'+str(BIAS(wrf_total,obs_total)))
        print( 'PEARSON TOTAL:',pearsonr(obs_total,wrf_total)[0])
        print( 'SPEARMAN TOTAL:',spearmanr(wrf_total,obs_total)[0])
        print( 'RMSE TOTAL', RMSE(wrf_total,obs_total))
        print( 'BIAS TOTAL', BIAS(wrf_total,obs_total))



if __name__ == "__main__":
    main()

