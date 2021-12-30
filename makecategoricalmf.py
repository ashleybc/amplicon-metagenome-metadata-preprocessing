import numpy as np
import pandas as pd
from collections import defaultdict

#functions
def contmapflproc(mf,of,od,basicInfoDict,omitList=[]): #mf is mapping file; of=oxygen classification file; omitList = samples to drop,od=non-root output directory

    toRem=[]
    for i in mf.index:
        if i in omitList or "blank" in i:
            toRem.append(i)

    #show the rows in metadata table that were not successfully sequenced
    mf.loc[toRem]

    #drop from metadata
    mf=mf.drop(toRem)

    #need to add new cols to cleaned mapping file and pre-fill with NaN for later
    newcols=["oxcat","redoxclinepos","Month","year"]
    
    for c in newcols:
        mf[c]=np.nan

    #loop through conditions. format of O2cond table is so that comparisons to mapfl data should be >= or <=
    for ind in of.index: #for each redox category in O2cond
        sulf=of["Sulfide"].loc[ind] #sulfide Y/N
        ubox=of["ubO2"].loc[ind] #>=DO upper boundary
        lbox=of["lbO2"].loc[ind]    #<=DO lower boundary
        print("lbox",lbox)
        oxname=ind

        if sulf=="Y": #if sample is euxinic
            locs=mf[mf["sulfide_uM"]!=0].index
            mf["oxcat"].loc[locs]=ind #set oxcat in mapfile
        
        else:
            locs=mf[(mf["DO_um"]>=lbox) & (mf["DO_um"]<=ubox)].index #find appropriate range of O2 conditions in mapfl
            mf["oxcat"].loc[locs]=ind #set oxcat in mapfl


    indList=list(mf.index)
    for k in list(basicInfoDict.keys()):
        print(k)
            
        checkList= [i for i in indList if k in i]
            
        if checkList:
            print(k)
            subdf=mf.loc[checkList]
            #print(subdf)

            above_locs=subdf[(subdf["Depth"]<basicInfoDict[k]["URB"])].index
            mf["redoxclinepos"].loc[above_locs]="above"

            URB_locs=subdf[(subdf["Depth"]==basicInfoDict[k]["URB"])].index
            mf["redoxclinepos"].loc[URB_locs]="upper boundary"
 
            inside_locs=subdf[(subdf["Depth"]>basicInfoDict[k]["URB"])&(subdf["Depth"]<basicInfoDict[k]["LRB"])].index
            mf["redoxclinepos"].loc[inside_locs]="inside"
                
            LRB_locs=subdf[(subdf["Depth"]==basicInfoDict[k]["LRB"])].index
            mf["redoxclinepos"].loc[LRB_locs]="lower boundary"
        
            below_locs=subdf[(subdf["Depth"]>basicInfoDict[k]["LRB"])].index
            mf["redoxclinepos"].loc[below_locs]="below"
                
    #fill basic info

    for i in mf.index:
        for k in list(basicInfoDict.keys()):
            if k in i:
                mf["Month"].loc[i]=basicInfoDict[k]["Month"]
                mf["year"].loc[i]=basicInfoDict[k]["year"]

    mf.to_csv(root+"/"+od+"/cleanedcontmapfl.csv",index=True)
    
    return mf


def catmapgen(mf,mfcc,mdbin,od): #mf=continuous mapping file, mfcc=skip columns; mdbin=metadata bins file; od=non-root output directory
    
    mfcopy=mf.copy()
    
    for col in mfcopy.columns:
        if col not in mfcc:
        
            col2 = col + "Cat"
            mfcopy[col2] = mfcopy[col]
            print(mfcopy[col])

            binseps = eval(mdbin.loc[col].bin_separators)
            
            for ind, be in enumerate(binseps):
            
                if ind == 0:
                    mfcopy[col2][mfcopy[col] < be] = str(ind)
                else:
                    mfcopy[col2][(mfcopy[col] < be) & (mfcopy[col] >= binseps[ind-1]) ] = str(ind)
                    
            mfcopy[col2][(mfcopy[col] >= binseps[-1])] = str(ind + 1)
        
            mfcopy[col2][np.isnan(mfcopy[col])] = "NAN"
            print(mfcopy[col2])

        elif col in mfcc:
        
            col2 = col + "Cat"
            mfcopy[col2] = mfcopy[col]
        
        else:
            print("ERROR")
            
    catCols=[col for col in mfcopy.columns if "Cat" in col]
    catmf=mfcopy[catCols]

    
    catmf.to_csv(root+"/"+od+"/cleanedcatmapfl.csv",index=True)
    
    return catmf

def watercol_runner():

    #metadata... need to also stick the del13POC data onto the original metadatafile
    mapfl=pd.read_csv("newmetadatafiles/AshFGL_mappingfile.tsv",delimiter="\t",index_col="#Sample ID")
    metadatabins=pd.read_csv("nonqiimeinput/metadatabinning.txt", sep="\t",index_col=0)
    
    #user inputs
    dropsamples=["AC-GL4-68","AC-GL4-68C","AC-GL4-55","AC-GL4-55C","AC-GL3-2-7-14","AC-GL3-0-2-14"]
    mfcatignorecol=["SizeFrac","Description","oxcat","redoxclinepos","Month","year","pos"]

    cleanedcontmapfl=contmapflproc(mapfl,o2cond,"newmetadatafiles",cruiseDict,omitList=dropsamples) #mf is mapping file; of=oxygen classification file; omitList = samples to drop,od=non-root output directory
    catmapfl=catmapgen(cleanedcontmapfl,mfcatignorecol,metadatabins,"newmetadatafiles")

def trap_runner():

    #metadata... need to also stick the del13POC data onto the original metadatafile
    mapfl=pd.read_csv("/Users/ashley/Documents/Research/Gordon\'sLab/FGL/amplicon/Qiime2_analyzed_FGL_traps/metadata.tsv",delimiter="\t",index_col="#Sample ID")
    metadatabins=pd.read_csv("/Users/ashley/Documents/Research/Gordon\'sLab/FGL/amplicon/Qiime2_analyzed_FGL_traps/trapmetadatabinning.txt", sep="\t",index_col=0)

    #user inputs
    dropsamples=["tr-blk","GEN-DONOR","MOCK-EVEN","MOCK-STAG"]
    mfcatignorecol=["oxcat","redoxclinepos","Month","year","pos"]
    
    cleanedcontmapfl=contmapflproc(mapfl,o2cond,"newmetadatafilestraps",cruiseDict,omitList=dropsamples) #mf is mapping file; of=oxygen classification file; omitList = samples to drop,od=non-root output directory
    catmapfl=catmapgen(cleanedcontmapfl,mfcatignorecol,metadatabins,"newmetadatafilestraps")

#body

#file import
root="/Users/ashley/Documents/Research/Gordon'sLab/FGL/amplicon/pythoncode"
o2cond=pd.read_csv("nonqiimeinput/O2classNew.csv",index_col="term")
#multiply infinities in anoxic/euxinic water to = -infinity
o2cond["lbO2"]=o2cond["lbO2"].apply(lambda x: x*-1 if x==np.inf else x)

cruiseDict=defaultdict(dict)

cruiseDict["GL1"]["Month"]="early July"
cruiseDict["GL2"]["Month"]="mid August"
cruiseDict["GL3"]["Month"]="early October"
cruiseDict["GL4"]["Month"]="late July"

for k in list(cruiseDict.keys()):
    if k !="GL4":
        cruiseDict[k]["year"]=2016
    else:
        cruiseDict[k]["year"]=2017
    
    if "July" in cruiseDict[k]["Month"]:
        cruiseDict[k]["LRB"]=20.0
    elif "October" in cruiseDict[k]["Month"]:
        cruiseDict[k]["LRB"]=21.0
    elif "August" in cruiseDict[k]["Month"]:
        cruiseDict[k]["LRB"]=15.0
        
    if "August" not in cruiseDict[k]["Month"]:
        cruiseDict[k]["URB"]=15.0
    else:
        cruiseDict[k]["URB"]=13.36
        
