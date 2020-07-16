# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:37:28 2020

@author: theodore
"""
import os.path
import os
import subprocess
from subprocess import call
import math
import itertools
from Bio.Blast.Applications import NcbiblastpCommandline

def checkKey(dict, key): 
      
    if key in dict.keys(): 
        #print("Present, ", end =" ") 
        #print("value =", dict[key]) 
        return True
    else: 
        #print("Not present") 
        return False




cwd = os.getcwd()
#Fastaname = '1a2y.fasta'
#PDBname = '1a2y.pdb'

Fastaname = '1a3l.fasta'
PDBname = '1a3l.pdb'

toolsPath = '/Tools/'
msmsPath = '/Tools/MSMS/'
writeFilePath = '/ProcessFile/'
processFileName = ""
PDBid = ""
blastShortPath='/Tools/Blast/Epitope_SV'
def extractPDBContent():
    global processFileName
    filename = cwd+writeFilePath+PDBname
    PDB = PDBname[0:-4]
    with open(filename) as f:
        contents = f.readlines()
    f.close()
    txt = ""
    for c in contents:
        if 'REMARK' in c:
            continue
        if "HELIX" in c:
            #print(c)
            txt+=(c)
        if "SHEET" in c:
            #print(c)
            txt+=(c)
        if "ATOM" in c:
            if "HOH" in c:
                continue
            if "PO4" in c:
                continue
            txt+=(c)
    #print(txt)
    processFileName = cwd+writeFilePath+PDB+"_new"+".pdb"
    with open(cwd+writeFilePath+PDB+"_new"+".pdb",'w') as f:
        f.write(txt)
    f.close()


def getPDBContentFromFile():
    global TotalAtomNum,TotalAANum,PDBallchain
    global ATOMdata,PDBFile
    TotalAtomNum = 0
    TotalAANum = 0
    ATOMdata = {}
    PDBFile = {}
    PDBallchain=''
    NowResNum = NowChain = ''
    #filename = cwd+writeFilePath+PDBname
    if os.path.exists(processFileName):
        print('file exist')
    else:
        print('No')
        return
    with open(processFileName) as f:
        contents = f.readlines()
    for Str in contents:
        Str = Str[:80]
        if len(Str) <= 1:
            continue
        token = Str[:6] 
        if token == 'ENDMDL':
            break
        if token == 'ATOM  ':
            #print(Str)
            AtomName = Str[12:16].strip(' ')
            AltLoc = Str[16:17]
            #print(AltLoc)
            if AltLoc == ' ' or AltLoc == 'A':
                ResName = Str[17:20].strip(' ')
            else:
                continue
            ChainID = Str[21:22].strip(' ')
            ResNum = Str[22:27].strip(' ')
            CoordX = (Str[30:38].strip(' '))
            CoordY = (Str[38:46].strip(' '))
            CoordZ = (Str[46:54].strip(' '))
            BFactor = (Str[60:66].strip(' '))
            if BFactor == '':
                BFactor = 0
            TotalAtomNum+=1
            if NowResNum!=ResNum or NowChain!=ChainID:
                TotalAANum+=1
                if NowChain!=ChainID:
                    PDBallchain += ChainID
                NowResNum = ResNum
                NowChain = ChainID
                if ResName == 'GLY':
                    ResName ='G'
                elif ResName == 'ALA':
                    ResName ='A'
                elif ResName == 'LEU':
                    ResName ='L'
                elif ResName == 'VAL':
                    ResName ='V'
                elif ResName == 'ILE':
                    ResName ='I'
                elif ResName == 'PHE':
                    ResName ='F'
                elif ResName == 'TRP':
                    ResName ='W'
                elif ResName == 'TYR':
                    ResName ='Y'
                elif ResName == 'ASP':
                    ResName ='D'
                elif ResName == 'HIS':
                    ResName ='H'
                elif ResName == 'ASN':
                    ResName ='N'
                elif ResName == 'GLU':
                    ResName ='E'
                elif ResName == 'LYS':
                    ResName ='K'
                elif ResName == 'GLN':
                    ResName ='Q'
                elif ResName == 'MET':
                    ResName ='M'
                elif ResName == 'ARG':
                    ResName ='R'
                elif ResName == 'SER':
                    ResName ='S'
                elif ResName == 'THR':
                    ResName ='T'
                elif ResName == 'CYS':
                    ResName ='C'
                elif ResName == 'PRO':
                    ResName ='P'
                else:
                    print('cant identify the residue{} {}'.format(ResNum,ResName))
                    ResName = 'X'
                if TotalAANum<10:
                   print(AtomName,AltLoc,ResName,ChainID,ResNum,BFactor)
                if ChainID not in PDBFile.keys():
                   PDBFile[ChainID] = {}
                if ResName not in PDBFile[ChainID].keys():
                   PDBFile[ChainID][ResNum] = {}
                PDBFile[ChainID][ResNum]['Res'] = ResName
                PDBFile[ChainID][ResNum]['BFactor'] = BFactor
                PDBFile[ChainID][ResNum]['AtomNum'] = 0
            PDBFile[ChainID][ResNum]['AtomNum'] += 1
            if 'Pos' not in PDBFile[ChainID][ResNum].keys():
                PDBFile[ChainID][ResNum]['Pos'] = {}
            #print(AtomName,AltLoc,ResName,ChainID,ResNum,PDBFile[ChainID][ResNum]['Pos'])
            if AtomName not in PDBFile[ChainID][ResNum]['Pos'].keys():
                #print(AtomName,ChainID,ResNum,PDBFile[ChainID][ResNum]['Pos'])
                PDBFile[ChainID][ResNum]['Pos'][AtomName] = {}
            PDBFile[ChainID][ResNum]['Pos'][AtomName]['X'] = CoordX
            PDBFile[ChainID][ResNum]['Pos'][AtomName]['Y'] = CoordY
            PDBFile[ChainID][ResNum]['Pos'][AtomName]['Z'] = CoordZ
    #print(PDBFile)
    f.close()
    ATOMdata = PDBFile
    #print(ATOMdata)  
    #print(TotalAtomNum)
    #print(PDBFile)



#print(TotalAtomNum,TotalAANum,ATOMdata)
def ExeMSMS():
    msmsOutFile = ['face','vert','area']
    msmsoutputfileExistingFlag = 1;
    for fileType in msmsOutFile:
        print(msmsPath+'MSMS_file/'+Fastaname+'.'+fileType)
        if os.path.exists(msmsPath+'MSMS_file/'+Fastaname+'.'+fileType) == False:
            msmsoutputfileExistingFlag = 0
    if msmsoutputfileExistingFlag == 1:
        print('File Exist'+msmsPath+'MSMS_file/'+Fastaname+'.'+fileType)
        return 0;
    #File Non-exist then exec MSMS
    PDB = PDBname[:-4]
    upload_path = writeFilePath
    print(PDB+'   '+upload_path)
    MSMS_Path = '/Tools/MSMS/'
    MSMS_Output_file_path = '/Tools/MSMS/MSMS_file/'
    XYZ_File_Path = '/Tools/MSMS/XYZ_file/'
    
    if os.path.exists(cwd+"{}{}.pdb".format(upload_path,PDB)) == False:
        print('error in EXEMSMS cant found {}{}.pdb'.format(upload_path,PDB))
        #return
    else:
        print('PDB exist')
    if os.path.exists(cwd+MSMS_Path+'pdb_to_xyzrn.py') == False:
        print('error in EXEMSMS cant found'+MSMS_Path+'pdb_to_xyzrn.py')
    else:
        print('PDB_TO_xyzrn exist')
    p = subprocess.Popen('python Tools/MSMS/pdb_to_xyzrn.py ProcessFile/{}.pdb > Tools/MSMS/XYZ_file/{}.xyzrn'.format(PDB+"_new",PDB),shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    #print(out,err)
    if os.path.exists(cwd+MSMS_Path+'XYZ_file/{}.xyzrn'.format(PDB)) == False:
        print('error in EXEMSMS cant found'+MSMS_Path+'XYZ_file/{}.xyzrn'.format(PDB))
    else:
        print('xyzrn file exist')
    
    #print("msms.exe -if {}{}.xyzrn -af {} -of {} -no_header".format(XYZ_File_Path[1:],PDB,PDB,PDB,PDB))
    p = subprocess.Popen("msms.exe -if {}{}.xyzrn -af {} -of {} -no_header".format(XYZ_File_Path[1:],PDB,PDB,PDB,PDB),shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
#extractPDBContent()
#getPDBContentFromFile()  
#ExeMSMS()  
MSMS = {}
NeiAATable = {} #記得放回去
def intergrateMSMSInfo():
    #記得要去改Vert Face 這三個檔案的位置
    Path_tmp=msmsPath
    PDB = PDBname[:-4]
    with open(PDB+".vert") as f:
        contents = f.readlines()
    MSMS['VERT'] = {}
    VectNum = 1
    for i in contents:
        if VectNum not in MSMS['VERT']:
            MSMS['VERT'][VectNum] = {}
        line = i.split(' ')
        line = list(filter(None, line))
        x = (line[0]);
        y = (line[1]);
        z = (line[2]);
        tmp = line[-1].split('_')
        MSMS['VERT'][VectNum]['x'] = x
        MSMS['VERT'][VectNum]['y'] = y
        MSMS['VERT'][VectNum]['z'] = z
        MSMS['VERT'][VectNum]['ResAndChain'] = tmp[2].strip('\n')
        VectNum+=1
    f.close()
    #print(tmp[2].strip('\n'))
    #print(VectNum)
    #print(MSMS)
    with open(PDB+'.face') as f:
        contents = f.readlines()
    MSMS['FACE'] = {}
    FaceNum = 1
    for i in contents:
        if FaceNum not in MSMS['FACE']:
            MSMS['FACE'][FaceNum] = {}
        line = i.split(' ')
        line = list(filter(None, line))
        a = int(line[0]);
        b = int(line[1]);
        c = int(line[2]);
        MSMS['FACE'][FaceNum]['a'] = a
        MSMS['FACE'][FaceNum]['b'] = b
        MSMS['FACE'][FaceNum]['c'] = c
        #print('face',a,b,c)
        FaceNum+=1
        #print(a,b,c)
        #print(MSMS['VERT'][a]['ResAndChain'])
        a = MSMS['VERT'][a]['ResAndChain']
        b = MSMS['VERT'][b]['ResAndChain']
        c = MSMS['VERT'][c]['ResAndChain']
        #print('Vert',a,b,c)
        if (a==b and b == c):
            continue
        if (a != b):
            '''
             Python 不支援直接塞的功能 所以先check是否存在 在判定是否為空
            '''
            if checkKey(NeiAATable,a) == False:
                NeiAATable[a] = {}
            if checkKey(NeiAATable[a],b) == False:
                NeiAATable[a][b] = 1
            if NeiAATable[a][b] >= 1:
                NeiAATable[a][b]+=1
            
            if checkKey(NeiAATable,b) == False:
                NeiAATable[b] = {}
            if checkKey(NeiAATable[b],a) == False:
                NeiAATable[b][a] = 1
            if NeiAATable[b][a] >= 1:
                NeiAATable[b][a]+=1
        if (b != c):
            '''
             Python 不支援直接塞的功能 所以先check是否存在 在判定是否為空
            '''
            if checkKey(NeiAATable,b) == False:
                NeiAATable[b] = {}
            if checkKey(NeiAATable[b],c) == False:
                NeiAATable[b][c] = 1
            if NeiAATable[b][c] >= 1:
                NeiAATable[b][c]+=1
            
            if checkKey(NeiAATable,c) == False:
                NeiAATable[c] = {}
            if checkKey(NeiAATable[c],b) == False:
                NeiAATable[c][b] = 1
            if NeiAATable[c][b] >= 1:
                NeiAATable[c][b]+=1
        if (a != c):
            '''
             Python 不支援直接塞的功能 所以先check是否存在 在判定是否為空
            '''
            if checkKey(NeiAATable,a) == False:
                NeiAATable[a] = {}
            if checkKey(NeiAATable[a],c) == False:
                NeiAATable[a][c] = 1
            if NeiAATable[a][c] >= 1:
                NeiAATable[a][c]+=1
            
            if checkKey(NeiAATable,c) == False:
                NeiAATable[c] = {}
            if checkKey(NeiAATable[c],a) == False:
                NeiAATable[c][a] = 1
            if NeiAATable[c][a] >= 1:
                NeiAATable[c][a]+=1
    f.close()
    #print(MSMS)
    for key,NeiAAs in NeiAATable.items():
        AAtmp = key.split(':')
        #print(AAtmp[1],AAtmp[0])
        if AAtmp[0] in ATOMdata[AAtmp[1]].keys():
            #print(ATOMdata[AAtmp[1]][AAtmp[0]])
            ATOMdata[AAtmp[1]][AAtmp[0]]['NeiAAID'] = ''
            for NeiAA,Num in NeiAAs.items():
                ATOMdata[AAtmp[1]][AAtmp[0]]['NeiAAID'] += (str(NeiAA)+',')
            ATOMdata[AAtmp[1]][AAtmp[0]]['NeiAAID'] = ATOMdata[AAtmp[1]][AAtmp[0]]['NeiAAID'][:-1]
            #for v in NeiAA
        #B 15
        #print(ATOMdata[AAtmp[1]][AAtmp[0]])
        #ATOM
    #print(NeiAATable)
    """
    取得Area的表面面積的資訊
    """
    with open(PDB+".area") as f:
        contents = f.readlines()
    #MSMS['VERT'] = {}
    for line in contents:
        if len(line) == 0:
            continue
        SAS = line[16:24]
        AreaTmp = line[24:].split('_')
        if len(AreaTmp)!=3:
            continue
        AATmp = AreaTmp[2].split(':')
        Loc = AATmp[0]
        Chain = AATmp[1].strip('\n')
        #print(SAS,AATmp,Loc,Chain)
        if Loc in ATOMdata[Chain].keys():
           #print('AAA',ATOMdata[Chain][Loc])
           if 'msmsArea' not in ATOMdata[Chain][Loc].keys():
               ATOMdata[Chain][Loc]['msmsArea'] = SAS
           #print('AAA',ATOMdata[Chain][Loc]['msmsArea'])
           else:
               ATOMdata[Chain][Loc]['msmsArea'] += SAS
           #print('AAA',ATOMdata[Chain][Loc]['msmsArea'])
           #print(ATOMdata[Chain][Loc])
           if AreaTmp[0] in ATOMdata[Chain][Loc]['Pos'].keys():
               ATOMdata[Chain][Loc]['Pos'][AreaTmp[0]]['area'] = SAS
        #if AAtmp[0] not in ATOMdata[AAtmp[1]].keys():

def getTwoResidueShortDisTable(ResAname,ResBname):
    dis = 999999;
    ResAname = ResAname.split(':')
    ResBname = ResBname.split(':')
    if ResAname[0] not in ATOMdata[ResAname[1]].keys():
        return dis
    if ResBname[0] not in ATOMdata[ResBname[1]].keys():
        return dis
    ResA =  ATOMdata[ResAname[1]][ResAname[0]]
    ResB =  ATOMdata[ResBname[1]][ResBname[0]]
    
    for AtomType, AtomAinfo in ResA['Pos'].items():
        #print(AtomAinfo['area'])
        if 'area' not in AtomAinfo.keys():
            continue
        if float(AtomAinfo['area']) == 0.0:
            continue
        if ('X' not in AtomAinfo.keys()) or ('Y' not in AtomAinfo.keys()) or ('Z' not in AtomAinfo.keys()):
            
            continue
        x1 = float(AtomAinfo['X'])
        y1 = float(AtomAinfo['Y'])
        z1 = float(AtomAinfo['Z'])
        for AtomType, AtomBinfo in ResB['Pos'].items():
            if 'area' not in AtomBinfo.keys():
                continue
            if float(AtomBinfo['area']) == 0.0:
                continue
            if ('X' not in AtomBinfo.keys()) or ('Y' not in AtomBinfo.keys()) or ('Z' not in AtomBinfo.keys()):
                continue
            x2 = float(AtomBinfo['X'])
            y2 = float(AtomBinfo['Y'])
            z2 = float(AtomBinfo['Z'])
            tmp = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
            #print(tmp)        
            if dis > tmp:
                dis = tmp
    dis = math.sqrt(dis)
    return dis
        #print(AtomType,AtomAinfo)
    
        #print(AtomType,AtomBinfo)
def deleteIsNOTneighborResidues(NeiResIDs,ShortDisTable):
    AllResNum = len(NeiResIDs)
    for no,ResId in enumerate(NeiResIDs):
        #print(no,ResId)
        IsNotNeiResNum = 0
        for otherResID,Dis in ShortDistTable[ResId].items():
            #print(otherResID,Dis)
            if Dis == 999999:
                #print('yes')
                IsNotNeiResNum +=1
        #print(IsNotNeiResNum,AllResNum)
        if IsNotNeiResNum == AllResNum-1:
            #print('yes')
            NeiResIDs.remove(ResId)
            NeiResIDs = sorted(NeiResIDs)
            #print(NeiResIDs)
        #print('--------')
    return        
    
def bulitShortDisTableforNeighborAAs(NeiResIDs):
    NeiNum = len(NeiResIDs)
    NeiAADisTable = {}
    for i in range(0,NeiNum):
        for j in range(i+1,NeiNum):
            dis = getTwoResidueShortDisTable(NeiResIDs[i],NeiResIDs[j])
            if NeiResIDs[i] not in NeiAADisTable.keys():
                NeiAADisTable[NeiResIDs[i]] = {}
            #if NeiResIDs[j] not in NeiAADisTable[i].keys():
            #    NeiAADisTable[NeiResIDs[i]][NeiResIDs[j]] = {}    
            if NeiResIDs[j] not in NeiAADisTable.keys():
                NeiAADisTable[NeiResIDs[j]] = {}
            #if NeiResIDs[i] not in NeiAADisTable[j].keys():
            #    NeiAADisTable[NeiResIDs[j]][NeiResIDs[i]] = {}
            
            NeiAADisTable[NeiResIDs[i]][NeiResIDs[j]] = dis
            NeiAADisTable[NeiResIDs[j]][NeiResIDs[i]] = dis
            #print(dis)
    #print(NeiAADisTable)        
    return NeiAADisTable
'''
def combination(Dict,Res):
    if (not Dict) == True:
        print('asdf')
        #Dict[] = {}
    else:
        for k,v in Dict.items():
'''           
def correctSpiralSeq():
    for Chain,Chain_info in ATOMdata.items():
        for AA_resID,AA_ResInfo in Chain_info.items():
            if 'NeiAAID' not in ATOMdata[Chain][AA_resID].keys():
                continue
            else:
                if ATOMdata[Chain][AA_resID]['NeiAAID'] is None:
                    continue
            #print(ATOMdata[Chain][AA_resID]['NeiAAID'])
            NeiResIDs = ATOMdata[Chain][AA_resID]['NeiAAID'].split(',')
            #print(NeiResIDs)
            #建造鄰近胺基酸的最段距離
            ShortDistTable = bulitShortDisTableforNeighborAAs(NeiResIDs)
            #print(NeiResIDs)
            """
            deleteIsNOTneighborResidues 不能用 先直接寫 待檢討
            """
            AllResNum = len(NeiResIDs)
            for no,ResId in enumerate(NeiResIDs):
                #print(no,ResId)
                IsNotNeiResNum = 0
                #print(ShortDistTable[ResId])
                for otherResID,Dis in ShortDistTable[ResId].items():
                    if Dis == 999999:
                        #print('yes')
                        IsNotNeiResNum +=1
                #print(IsNotNeiResNum,AllResNum)
                if IsNotNeiResNum == AllResNum-1:
                     NeiResIDs.remove(ResId)
                     NeiResIDs = sorted(NeiResIDs)
                     #print(NeiResIDs)
                #print('----------')
            """
            ====================deleteIsNOTneighborResidues結束=============
            """
            NeiNum = len(NeiResIDs)
            if NeiNum <=7:
                First = NeiResIDs[0]
                NeiResIDs.remove(First)   
                if NeiNum == 1:
                    continue
                #print(NeiResIDs)
                """
                 學長的combination不能直接改成python坂 先用itertools頂住
                """
                #print(NeiResIDs)
                NeiResIDs = list(itertools.permutations(NeiResIDs))
                
                #print(NeiResIDs)
                TSP_Dis = 999999;
                TSP_no = -1;
                for no,List in enumerate(NeiResIDs):
                    DisTmp = 0
                    List = [First]+list(List)+[First]
                    for j in range(NeiNum):
                        a = List[j]
                        b = List[j+1]
                        if a not in ShortDistTable.keys():
                            DisTmp = 999999
                            break
                        if b not in ShortDistTable[a].keys():
                            DisTmp = 999999
                            break
                        DisTmp += ShortDistTable[a][b]
                    if DisTmp < TSP_Dis:
                        TSP_Dis = DisTmp
                        TSP_no = no
                ATOMdata[Chain][AA_resID]['NeiAAID'] = [First]+list(NeiResIDs[TSP_no])
            else: #如果>7就不用暴力解
                #print("out of 7")
                Seqtmp = str(NeiResIDs[0])
                #print(Seqtmp)
                for j in range(NeiNum-2):
                    selectAA = j+1
                    DisTmp = 999999
                    #print(selectAA)
                    for k in range(j+1,NeiNum):
                        #print(ShortDistTable[NeiResIDs[j]][NeiResIDs[k]])
                        if DisTmp > ShortDistTable[NeiResIDs[j]][NeiResIDs[k]]:
                            DisTmp = ShortDistTable[NeiResIDs[j]][NeiResIDs[k]]
                            selectAA = k
                    tmp = NeiResIDs[selectAA]
                    NeiResIDs[selectAA ] = NeiResIDs[j+1]
                    NeiResIDs[j+1] = tmp
                    Seqtmp += (","+NeiResIDs[j+1])
                Seqtmp = Seqtmp.split(',')
                #print(Seqtmp)
                ATOMdata[Chain][AA_resID]['NeiAAID'] = Seqtmp.append(NeiResIDs[NeiNum-1])
                
                #ATOMdata[Chain][AA_resID]['NeiAAID'] = Seqtmp
                    #print('======='+str(DisTmp)+" "+str(TSP_Dis)+" "+str(TSP_no));
                #print(ATOMdata[Chain][AA_resID]['NeiAAID'])
                #print([First]+list(NeiResIDs[TSP_no]))
                
            #deleteIsNOTneighborResidues(NeiResIDs,ShortDistTable)
            #print(ShortDistTable)

def convertRes(ResName):
    if ResName == 'GLY':
        ResName ='G'
    elif ResName == 'ALA':
        ResName ='A'
    elif ResName == 'LEU':
        ResName ='L'
    elif ResName == 'VAL':
        ResName ='V'
    elif ResName == 'ILE':
        ResName ='I'
    elif ResName == 'PHE':
        ResName ='F'
    elif ResName == 'TRP':
        ResName ='W'
    elif ResName == 'TYR':
        ResName ='Y'
    elif ResName == 'ASP':
        ResName ='D'
    elif ResName == 'HIS':
        ResName ='H'
    elif ResName == 'ASN':
        ResName ='N'
    elif ResName == 'GLU':
        ResName ='E'
    elif ResName == 'LYS':
        ResName ='K'
    elif ResName == 'GLN':
        ResName ='Q'
    elif ResName == 'MET':
        ResName ='M'
    elif ResName == 'ARG':
        ResName ='R'
    elif ResName == 'SER':
        ResName ='S'
    elif ResName == 'THR':
        ResName ='T'
    elif ResName == 'CYS':
        ResName ='C'
    elif ResName == 'PRO':
        ResName ='P'
    else:
        print('cant identify the residue{}'.format(ResName))
        ResName = 'X'
    return ResName

from Bio.PDB.PDBParser import PDBParser
ID_table = dict()
def getPDBID():
    filename = cwd+writeFilePath+PDBname
    parser = PDBParser()
    structure = parser.get_structure('tmp',filename)
    model = structure[0]
    for chain in model:
        print(chain.id)
        ID_table[chain.id] = {}
        for i in chain:
            res = i.get_resname()
            ID = i.id
            if res == 'HOH' or res == 'POS':
                continue
            ID_table[chain.id][str(ID[1])] = convertRes(res)

getPDBID()    
#A 107
#B 116
#C 129
def getResidueType(Res):
    Res = Res.split(':')
    #print(Res)
    return ID_table[Res[1]][Res[0]]


def getNeighborResidueTypeList(NeiAAID):
    if type(NeiAAID) == str:
        NeiAAID = NeiAAID.split(',')
    #print(NeiAAID)
    AAlist = ""
    for res in NeiAAID:
        #print('res',res)
        try:
            AAlist += getResidueType(res)
        except:
            #
            pass
    return AAlist
    

#將表面胺基酸之各個螺旋特徵寫成fasta檔案
def writeSpiralVectorFastafileforBlastpShort():
    global PDBid
    PDBid = PDBname[0:-4]
    FastaName = PDBid+'_2allsurfaceresiduesfasta.fa'
    TxT = ""
    for Chain,Chain_info in ATOMdata.items():
        for AA_resID,AA_resInfo in Chain_info.items():
            if 'NeiAAID' not in AA_resInfo.keys():
                continue
            if AA_resInfo['NeiAAID'] is None:
                continue
            if AA_resInfo['msmsArea'] == 0:
                continue
            if AA_resInfo['AtomNum'] is None:
                continue
            print(AA_resInfo['NeiAAID'])
            
            #$ResidueOutputInfo = "{$PDBid}_{$Chain}_{$AA_resID}_{$AA_resInfo['Res']}";
            ResidueOutputInfo = "{}_{}_{}_{}".format(PDBid,Chain,AA_resID,AA_resInfo['Res'])
            #print(ResidueOutputInfo)
            TxT += (">"+ResidueOutputInfo+"_clockwise\n")
            SpiralVectorAAList = getNeighborResidueTypeList(AA_resInfo['NeiAAID'])
            print(SpiralVectorAAList)
            TxT += (SpiralVectorAAList+'\n')
            TxT += (">"+ResidueOutputInfo+"_anticlockwise\n")
            AntiSpiralVectorAAList = SpiralVectorAAList[::-1]
            #print(SpiralVectorAAList)
            TxT += (AntiSpiralVectorAAList+'\n')
    with open(FastaName,'w') as f:
        f.write(TxT)
    f.close()

def Blast():
    cwd = os.getcwd()
    FastaName = PDBid+'_2allsurfaceresiduesfasta.fa'
    print(cwd)
    print(FastaName)
    my_blast_exe = r"C:\Program Files\NCBI\blast-2.7.1+\bin\blastp.exe"
    BlastPath = cwd+blastShortPath
    OutputFile = 'blastshort.'+PDBname+'.xml'
    blastp_cline = NcbiblastpCommandline(cmd=my_blast_exe, query=FastaName, db=BlastPath,task = "blastp-short", evalue=1000,outfmt=5, out=OutputFile)
    print(blastp_cline)
    print(blastp_cline())

def extract_HitScore(Hsp):
    txt = ""
    for h in Hsp:
        if "bit" in str(h):
            txt = h.text
            break
    return txt
def extract_queryFrom(Hsp):
    txt = ""
    for h in Hsp:
        if "query-from" in str(h):
            txt = h.text
            break
    return txt
def extract_queryto(Hsp):
    txt = ""
    for h in Hsp:
        if "query-to" in str(h):
            txt = h.text
            break
    return txt
def extract_hitFrom(Hsp):
    txt = ""
    for h in Hsp:
        if "hit-from" in str(h):
            txt = h.text
            break
    return txt
def extract_hitto(Hsp):
    txt = ""
    for h in Hsp:
        if "hit-to" in str(h):
            txt = h.text
            break
    return txt
def extract_queryframe(Hsp):
    txt = ""
    for h in Hsp:
        if "query-frame" in str(h):
            txt = h.text
            break
    return txt

def extract_hitframe(Hsp):
    txt = ""
    for h in Hsp:
        if "hit-frame" in str(h):
            txt = h.text
            break
    return txt

def extract_alignlen(Hsp):
    txt = ""
    for h in Hsp:
        if "align-len" in str(h):
            txt = h.text
            break
    return txt
from bs4 import BeautifulSoup
def readXML():
    File = 'blastshort.'+PDBname+'.xml'
    soup = BeautifulSoup(open(File), "html.parser")
    Hits = soup.find_all('hit')
    Hits = soup.find_all('hit')
    print(len(Hits))
    Hits_array = {}
    cnt = 0
    for H in Hits:
            Hit_id = H.hit_id.text
            Hit_def = H.hit_def.text
            Hit_accession = H.hit_accession.text
            Hit_len = H.hit_len.text
            Hsp_array = {}
            cnt_hsp = 0
            for hsp in H.find_all('hsp'):
                hsp_num = hsp.hsp_num.text
                hsp_bit_score = extract_HitScore(hsp)
                hsp_score = hsp.hsp_score.text
                hsp_evalue = hsp.hsp_evalue.text
                hsp_query_from = extract_queryFrom(hsp)
                hsp_query_to = extract_queryto(hsp)
                hsp_hit_from = extract_hitFrom(hsp)
                hsp_hit_to = extract_hitto(hsp)
                hsp_query_frame = extract_queryframe(hsp)
                hsp_hit_frame = extract_hitframe(hsp)
                hsp_identity = hsp.hsp_identity.text
                hsp_positive = hsp.hsp_positive.text
                hsp_gaps = hsp.hsp_gaps.text
                try:
                    hsp_align_len = hsp.extract_alignlen(hsp)
                except:
                    hsp_align_len = 0
                hsp_qseq = hsp.hsp_qseq.text
                hsp_hseq = hsp.hsp_hseq.text
                hsp_midline = hsp.hsp_midline.text
                Hsp_array[cnt_hsp] = {"hsp_num":float(hsp_num),
                                   "hsp_bit_score":float(hsp_bit_score),
                                   "hsp_score":float(hsp_score),
                                   "hsp_evalue":float(hsp_evalue),
                                   "hsp_query_from":int(hsp_query_from),
                                   "hsp_query_to":int(hsp_query_to),
                                   "hsp_hit_from":int(hsp_hit_from),
                                   "hsp_hit_to":int(hsp_hit_to),
                                   "hsp_query_frame":int(hsp_query_frame),
                                   "hsp_hit_frame":int(hsp_hit_frame),
                                   "hsp_identity":int(hsp_identity),
                                   "hsp_positive":int(hsp_positive),
                                   "hsp_gaps":int(hsp_gaps),
                                   "hsp_align_len":int(hsp_align_len),
                                   "hsp_qseq":str(hsp_qseq),
                                   "hsp_hseq":str(hsp_hseq),
                                   "hsp_midline":str(hsp_midline),
                                  }
                cnt_hsp += 1
            Hits_array[cnt] = {
                              "Hit_id":str(Hit_id),
                              "Hit_def":str(Hit_def),
                              "Hit_accession":str(Hit_accession),
                              "Hit_len":str(Hit_len),
                              "Hit_hsps":Hsp_array,
                              }
            cnt+=1
        
def readXML2():
    File = 'blastshort.'+PDBname+'.xml'
    soup = BeautifulSoup(open(File), "html.parser")
    Iter = soup.find_all('iteration')
    Iter_Hits = (soup.find_all('iteration_hits'))            
    Iter_stats = (soup.find_all('iteration_stat'))
    Iter_message = (soup.find_all('iteration_message'))
    XML_data = dict()
    #print(len(Iter_Hits),len(Iter_stats),len(Iter_message))
    for idx,_ in enumerate(Iter_Hits):
        Hits = Iter_Hits[idx].find_all('hit')
        #print(len(Hits))
        Hits_array = {}
        cnt = 0
        for H in Hits:
                Hit_id = H.hit_id.text
                Hit_def = H.hit_def.text
                Hit_accession = H.hit_accession.text
                Hit_len = H.hit_len.text
                Hsp_array = {}
                cnt_hsp = 0
                for hsp in H.find_all('hsp'):
                    hsp_num = hsp.hsp_num.text
                    hsp_bit_score = extract_HitScore(hsp)
                    hsp_score = hsp.hsp_score.text
                    hsp_evalue = hsp.hsp_evalue.text
                    hsp_query_from = extract_queryFrom(hsp)
                    hsp_query_to = extract_queryto(hsp)
                    hsp_hit_from = extract_hitFrom(hsp)
                    hsp_hit_to = extract_hitto(hsp)
                    hsp_query_frame = extract_queryframe(hsp)
                    hsp_hit_frame = extract_hitframe(hsp)
                    hsp_identity = hsp.hsp_identity.text
                    hsp_positive = hsp.hsp_positive.text
                    hsp_gaps = hsp.hsp_gaps.text
                    try:
                        hsp_align_len = hsp.extract_alignlen(hsp)
                    except:
                        hsp_align_len = 0
                    hsp_qseq = hsp.hsp_qseq.text
                    hsp_hseq = hsp.hsp_hseq.text
                    hsp_midline = hsp.hsp_midline.text
                    Hsp_array[cnt_hsp] = {"hsp_num":float(hsp_num),
                                       "hsp_bit_score":float(hsp_bit_score),
                                       "hsp_score":float(hsp_score),
                                       "hsp_evalue":float(hsp_evalue),
                                       "hsp_query_from":int(hsp_query_from),
                                       "hsp_query_to":int(hsp_query_to),
                                       "hsp_hit_from":int(hsp_hit_from),
                                       "hsp_hit_to":int(hsp_hit_to),
                                       "hsp_query_frame":int(hsp_query_frame),
                                       "hsp_hit_frame":int(hsp_hit_frame),
                                       "hsp_identity":int(hsp_identity),
                                       "hsp_positive":int(hsp_positive),
                                       "hsp_gaps":int(hsp_gaps),
                                       "hsp_align_len":int(hsp_align_len),
                                       "hsp_qseq":str(hsp_qseq),
                                       "hsp_hseq":str(hsp_hseq),
                                       "hsp_midline":str(hsp_midline),
                                      }
                    cnt_hsp += 1
                Hits_array[cnt] = {
                                  "Hit_id":str(Hit_id),
                                  "Hit_def":str(Hit_def),
                                  "Hit_accession":str(Hit_accession),
                                  "Hit_len":str(Hit_len),
                                  "Hit_hsps":Hsp_array,
                                  }
                cnt+=1
        Stats = Iter_stats[idx].find('statistics')
        stats_array = dict()
        #print(Stats)
        for itm in Stats:
            if "-num" in str(itm):
                stats_array["Statistics_db-num"] = itm.text
            if "db-len" in str(itm):
                stats_array["Statistics_db-len"] = itm.text
            if "hsp-len" in str(itm):
                stats_array["Statistics_hsp-len"] = itm.text  
            if "-space" in str(itm):
                stats_array["Statistics_eff-space"] = itm.text     
        stats_array["Statistics_kappa"] = Stats.find('statistics_kappa').text
        stats_array["Statistics_lambda"] = Stats.find('statistics_lambda').text
        stats_array["Statistics_entropy"] = Stats.find('Statistics_entropy')
        It_message = Iter[idx].find('iteration_message')
        if It_message is None:
            iter_message = ""
        else:
            iter_message = It_message.text
        XML_data[idx] = {}
        for itr_itm in Iter[idx]:
            if "query-ID" in str(itr_itm):
                XML_data[idx]["ID"] = itr_itm.text
            if "query-def" in str(itr_itm):
                XML_data[idx]["query"] = itr_itm.text
            if "query-len" in str(itr_itm):
                XML_data[idx]["length"] = itr_itm.text
        XML_data[idx]["Hits"] = Hits_array
        XML_data[idx]["Statistics"] = stats_array
        XML_data[idx]["Iteration_message"] = iter_message

        #XML_data[idx]["ID"] = Iter[idx].find("")
        #print('iter_message',iter_message)
        #stats_array["Statistics_db-num"]
        #print(num,db_len,hsp_len,space,kappa,Lambda,entropy)
        
        #print(len(Stats))
        
    return XML_data    
        #print(len(Hits))
#readXML() 
def MappingEpitopeSpiralVectorPatch():
    for idx,SingleResidueInfo in BlastShortdata.items():
        if "No hits found" in SingleResidueInfo['Iteration_message']:
            continue
        Res = SingleResidueInfo['query'].split('_')
        Obj_chain = Res[1]
        Obj_pdbid = Res[2]
        #print(Obj_chain,Obj_pdbid,SingleResidueInfo['Hits'])
        for i,SV_Res_Hit in SingleResidueInfo['Hits'].items():
            HitSVResInfo = SV_Res_Hit['Hit_id']
            Evalue = SV_Res_Hit['Hit_hsps'][0]['hsp_evalue']
            print(HitSVResInfo,Evalue)
            print('-------')
            if 'Predict_Epitope' not in ATOMdata[Obj_chain][Obj_pdbid].keys():
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope'] = {}
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'] = {}
            if HitSVResInfo not in ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'].keys():
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] = 99999
            if ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] == 99999:
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] = Evalue
            else:
                if ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] > Evalue:
                    ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] = Evalue
            '''
            if 'Predict_Epitope' not in ATOMdata[Obj_chain][Obj_pdbid].keys():
                print('no')
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope'] = {}
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'] = {}
                ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] = Evalue
            else:
                print('yes')
                if HitSVResInfo not in ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'].keys:
                    print('no2')
                    ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'] = {}
                    
                if ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] > Evalue:
                    print('yes2')
                    ATOMdata[Obj_chain][Obj_pdbid]['Predict_Epitope']['similarSVP'][HitSVResInfo] = Evalue
            '''
        #print(Evalue)
    #for SV_Res_Hit in SingleResidueInfo['Hits']:
    #    HitSVResInfo = SV_Res_Hit['Hit_id']
#MappingEpitopeSpiralVectorPatch()
def deletesingleMappingResiduebySVP(killID):
    killID = killID.split(':')
    del ATOMdata[killID[1]][killID[0]]["Predict_Epitope"]["similarSVP"][killID[2]]
    Pred_Epi_SVP = ATOMdata[killID[1]][killID[0]]["Predict_Epitope"]["similarSVP"]
    if len(Pred_Epi_SVP) == 0:
        del ATOMdata[killID[1]][killID[0]]["Predict_Epitope"]["similarSVP"]
def getHitProteinListforSVP():
    PDBLists = dict()
    existList = dict()
    for chain,Chain_Info in ATOMdata.items():
        for AA_resID,AA_resInfo in Chain_Info.items():
            NowAA = AA_resID+':'+chain
            if 'Predict_Epitope' in AA_resInfo.keys():
                for hitPDBandAA,e_value in AA_resInfo['Predict_Epitope']['similarSVP'].items():
                    hitPDB = hitPDBandAA.split('_')
                    hitPDB = hitPDB[0]+'_'+hitPDB[1]
                    if hitPDB in PDBLists.keys():
                        AAnum = len(PDBLists[hitPDB])
                    else:
                        AAnum = 0
                    if hitPDB in existList.keys():
                        if NowAA in existList[hitPDB].keys():
                            mappingNum = existList[hitPDB][NowAA]
                            if PDBLists[hitPDB][mappingNum]["Evalue"] > e_value:
                                deletesingleMappingResiduebySVP(PDBLists[hitPDB][mappingNum]["Mapping"])
                                PDBLists[hitPDB][mappingNum]["Mapping"] = NowAA+":"+hitPDBandAA
                                PDBLists[hitPDB][mappingNum]["Evalue"] = e_value
                            else:
                                deletesingleMappingResiduebySVP(PDBLists[hitPDB][mappingNum]["Mapping"])
                        else:
                            if hitPDB not in existList.keys():
                                existList[hitPDB] = {}
                                existList[hitPDB][NowAA] = AAnum
                            if hitPDB not in PDBLists.keys():
                                PDBLists[hitPDB] = {}
                                PDBLists[hitPDB][AAnum] = {}
                                PDBLists[hitPDB][AAnum]["Mapping"] = NowAA
                                PDBLists[hitPDB][AAnum]["Evalue"] = e_value
                                

                                
            
def deletefewResGroupandExpendCombineNeighborResEpitope():
    #print('i dont want to do it')
    getHitProteinListforSVP()

extractPDBContent()
ExeMSMS()
getPDBContentFromFile()
intergrateMSMSInfo()
correctSpiralSeq()
#print(ATOMdata)
writeSpiralVectorFastafileforBlastpShort()

'''
extractPDBContent()
ExeMSMS()
getPDBContentFromFile()
intergrateMSMSInfo()
correctSpiralSeq()
writeSpiralVectorFastafileforBlastpShort()
Blast()
BlastShortdata = readXML2()
MappingEpitopeSpiralVectorPatch()
'''

#第三種預測方式(Feature)
# Solid Angle
'''
def ExeSA():
    if os.path.exists(cwd+writeFilePath+PDBid+".pdb") == False:
        print("Not Found")
        return 
    #p = subprocess.Popen('python Tools/MSMS/pdb_to_xyzrn.py ProcessFile/{}.pdb > Tools/MSMS/XYZ_file/{}.xyzrn'.format(PDB+"_new",PDB),shell=True, stdout=subprocess.PIPE)
'''


'''
def calculateAndRead_SAFeature():
    ExeSA()
    ReadSA()
'''

#ProSA
ProsaWinSize = 0
PDBid = '1a2y'
def ExeProsa2003(chainId,WinSize = 0):
    if chainId not in PDBallchain:
        print("Chain not in PDB file")
        return 0
    moveEneFile = "ENG_"+PDBid+"_"+chainId+"_"+str(WinSize)+".pdb";
    moveTxtFile = "ENG_"+PDBid+"_"+chainId+"_list.txt";
    PdbID = PDBid+'_'+chainId
    filename = "instruction_"+PdbID+"_"+str(WinSize)+".cmd"
    #filename = "C:\\ProSa2003\\bin\\"+"instruction_"+PdbID+"_"+str(WinSize)+".cmd"
    with open(filename,'w') as file:
        tmp = "read pdb {}.pdb,{} {}\n".format(PDBid,chainId,PdbID)
        file.write(tmp)
        tmp = "print residue mapping {} {}\n".format(PdbID,moveTxtFile)
        file.write(tmp)
        file.write("analyse energy *\n")
        file.write("winsize * {}\n".format(WinSize))
        tmp = "print energy {} {}\n".format(PdbID,moveEneFile)
        file.write(tmp)
        file.write("quit")
        file.close()
    prosaProgramPath="C:/ProSa2003/bin/prosa2003.bat"
    if os.path.exists(prosaProgramPath) == False:
        print('ProSa Not Exists')
        return 0
    if os.path.exists(filename) == False:
        print('filename Not Exists')
        return 0
    cmd  = "{} {}".format(prosaProgramPath,filename)
    print('-------------')
    print(cmd)
    print('-------------')
    p = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    print(out,err)
#ExeProsa2003("B")
def readMappingTable(PDBid,chainID):
    MapTable = dict()
    print("ENG_"+PDBid+"_"+chainID+"_list.txt")
    try:
        with open("ENG_"+PDBid+"_"+chainID+"_list.txt") as f:
            contents = f.readlines()
            for c in contents:
                if 'Object' in c:
                    continue
                if 'Prosa' in c:
                    continue
                strings = c.split(' ')
                strings[:] = [x for x in strings if x]
                #print(strings)
                MapTable[strings[0]] = strings[1]
        f.close()
        #print(MapTable)
        return MapTable
    except:
        print('can not open ENG_file')
        return
#readMappingTable(PDBid,'B')
def ReadProsa2003(AAchain,WinSize = 0):
    MapTable = readMappingTable(PDBid,AAchain)
    try:
        filePath = "ENG_"+PDBid+"_"+AAchain+"_"+str(WinSize)+".pdb.ana"
        with open(filePath) as file:
            contents = file.readlines()
            for c in contents:
                if 'Protein' in c:
                    continue
                if 'Window' in c:
                    continue
                strings = c.split(' ')
                strings[:] = [x for x in strings if x]
                ProsaID = strings[0]
                Eng_1 = strings[1]
                Eng_2 = strings[2]
                Eng_3 = strings[3].strip('\n')
                #
                AA_id = MapTable[ProsaID]
                print(ProsaID,AA_id,Eng_1,Eng_2,Eng_3)
                if AAchain not in ATOMdata.keys():
                    print('error in function AAchain not in ATOMdata')
                    continue
                elif AA_id not in ATOMdata[AAchain] :
                    print('error in function AA_id not in ATOMdata')
                    continue
                else:
                    print('A')
                    #ATOMdata[AAchain][AA_id] = {}
                    ATOMdata[AAchain][AA_id]['FeatureCollection'] = {}
                    ATOMdata[AAchain][AA_id]['FeatureCollection']['Prosa_ENG'] = Eng_3
                    #print(ATOMdata)
                
    except:
        print('cant open ana_file')
        return
#ReadProsa2003("B")

def calculateAndRead_EnergyFeature():
    WinSize = ProsaWinSize
    ChainList = [i for i in PDBallchain]
    for chain in ChainList:
        ExeProsa2003(chain,WinSize)
        ReadProsa2003(chain,WinSize)
'''
calculateAndRead_EnergyFeature      
'''
