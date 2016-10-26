######################################
##   FMD Serotyping Analysis        ##
######################################


#  1. Create New Result Table w/ (Position, Serotype, A, T, C, G, Probe ID)w/ unique name (input)FG

def Create_Results_Table(ResultTableName, OutputFileName, START):
    import MySQLdb

    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    print "Initializing Database..."
    if START == "TRUE":
        
        dropline00 = ("DROP TABLE IF EXISTS FMD_GPR")
        cursor.execute(dropline00)
        
        executeline00 = ("""
            CREATE TABLE FMD_GPR
            (
            Block_ DOUBLE,
            Column_ DOUBLE,
            Row DOUBLE,
            Name_ VARCHAR(250),
            ID_ VARCHAR(250),
            X_ DOUBLE,
            Y_ DOUBLE,
            Dia_ DOUBLE,
            F532_Median DOUBLE,
            F532_Mean DOUBLE,
            F532_SD DOUBLE,
            F532_CV DOUBLE,
            B532 DOUBLE,
            B532_Median DOUBLE,
            B532_Mean DOUBLE,
            B532_SD DOUBLE,
            B532_CV DOUBLE,
            PCT_B532_1SD DOUBLE,
            PCT_B532_2SD DOUBLE,
            F532pct_Sat DOUBLE,
            RatioOfMeans1_532 DOUBLE,
            MedianOfRatios1_532 DOUBLE,
            MeanOfRatios1_532 DOUBLE,
            RatioSD1_532 DOUBLE,
            RgnRatio1_532 DOUBLE,
            RgnR21_532 DOUBLE,
            F_Pixels DOUBLE,
            B_Pixels DOUBLE,
            Circularity DOUBLE,
            SumOfMedian1_532 DOUBLE,
            SumOfMeans1_532 DOUBLE,
            LogRatio DOUBLE,
            F532_Median_B532 DOUBLE,
            F532_Mean_B532 DOUBLE,
            F532_Total_Intensity DOUBLE,
            SNR_532 DOUBLE,
            Flags DOUBLE,
            Normalize DOUBLE,
            Autoflag DOUBLE,
            RefNumber DOUBLE,
            ControlType VARCHAR(250),
            GeneName VARCHAR(250),
            TopHit VARCHAR(250),
            Description VARCHAR(250),
            Sequence VARCHAR(200))""")
        
        cursor.execute(executeline00)
        
        dropline0 =("DROP TABLE IF EXISTS FMD_%s") % (ResultTableName)
        cursor.execute(dropline0)  
        executeline0 = ("""
                CREATE TABLE FMD_%s
                (
                ID_Num VARCHAR(10),
                ID_Name VARCHAR(250),
                Mean532 VARCHAR(10),
                Probe_Seq VARCHAR(200))
                """) % (ResultTableName)
     
        cursor.execute (executeline0) 
        dropline8 =("DROP TABLE IF EXISTS WindowTemplate")
        cursor.execute(dropline8)
    
        executeline8 = ("""CREATE TABLE WindowTemplate
                        (
                        StartPos INTEGER,
                        EndPos INTEGER,
                        Title VARCHAR(255),
                        Sequence VARCHAR(1500))                        
                        """)
        cursor.execute (executeline8)

    dropline1 =("DROP TABLE IF EXISTS RESULTS_%s") % (ResultTableName)
    cursor.execute(dropline1)  
    executeline1 = ("""
            CREATE TABLE RESULTS_%s
            (
            Position VARCHAR(5),
            Probe_ID VARCHAR(10),
            Nuc_Score VARCHAR(10),
            Nucleotide VARCHAR(4))
            """) % (ResultTableName)
 
    cursor.execute (executeline1) 

    dropline2 =("DROP TABLE IF EXISTS Prelim_Scoring_%s") % (ResultTableName)
    cursor.execute(dropline2)
    executeline2 = ("""
            CREATE TABLE Prelim_Scoring_%s
            (
            Count VARCHAR(5),
            GI VARCHAR(15),
            Probe_ID VARCHAR(125),
            Probe_Length VARCHAR(5),
            Identities VARCHAR(5))
            """) % (ResultTableName)

    cursor.execute (executeline2)

    dropline3 =("DROP TABLE IF EXISTS Scored_Probes_%s") % (ResultTableName)
    cursor.execute(dropline3)
    executeline3 = ("""
            CREATE TABLE Scored_Probes_%s
            (
            GI VARCHAR(15),
            Probe_ID VARCHAR(125),
            Ratio VARCHAR(10))
            """) % (ResultTableName)

    cursor.execute (executeline3)

    dropline4 =("DROP TABLE IF EXISTS RESULTS_preScore_%s") % (ResultTableName)
    cursor.execute(dropline4)  
    executeline4 = ("""
            CREATE TABLE RESULTS_preScore_%s
            (
            Position INTEGER,
            SumOfScore INTEGER,
            Nucleotide VARCHAR(4))
            """) % (ResultTableName)
    
    cursor.execute (executeline4)

    dropline5 =("DROP TABLE IF EXISTS RESULTS_preScore2_%s") % (ResultTableName)
    cursor.execute(dropline5)  
    executeline5 = ("""
            CREATE TABLE RESULTS_preScore2_%s
            (
            Position INTEGER,
            MaxOfScore INTEGER)
            """) % (ResultTableName)
    
    cursor.execute (executeline5)

    dropline6 =("DROP TABLE IF EXISTS RESULTS_SCORED_%s") % (ResultTableName)
    cursor.execute(dropline6)  
    executeline6 = ("""
            CREATE TABLE RESULTS_SCORED_%s
            (
            Position INTEGER,
            Nucleotide VARCHAR(4),
            MFI INTEGER)
            """) % (ResultTableName)
    
    cursor.execute (executeline6)

    dropline7 =("DROP TABLE IF EXISTS RESULTS_PREconsensus_%s")%(ResultTableName)
    cursor.execute(dropline7)
    executeline7 = ("""
            CREATE TABLE RESULTS_PREconsensus_%s
            (
            Position INTEGER,
            MaxMFI INTEGER)            
            """)% (ResultTableName) #<< Added table for building final Consensus

    cursor.execute (executeline7)
    


    FileOUT_Name = """/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/%s.FASTA""" % (OutputFileName)
    
    Result_File = open(FileOUT_Name, "w")
    Result_File.close

def FASTACMDfile(GI,dbfile,outputfile):
    import os
    import sys
    
    FASTACMD_CMDLINE = "%sblastdbcmd -db %s -dbtype nucl -entry %s -target_only -out %s"
    status = os.system(FASTACMD_CMDLINE % (UtilDir,dbfile, GI, outputfile))

def FormatDB(Inputfile, proteinDB_TF, DBtitle):
    import os
    import sys
    try:
        print "Building BLAST database..."
        FORMATDB_CMDLINE = "/users/rwbarrettemac/bioinformatics/NCBIlegacy/bin/formatdb -i %s -p %s -o T -n %s"
        status = os.system(FORMATDB_CMDLINE % (UtilDir,Inputfile, proteinDB_TF, DBtitle))
    except:
        print "DB build failed..."

def MakeDB(Inputfile, dbType, DBtitle):
    import os
    import sys
    
    print "Building BLAST database..."
    MAKEDB_CMDLINE = "%smakeblastdb -in %s -out %s -parse_seqids -dbtype %s"
    status = os.system(MAKEDB_CMDLINE % (UtilDir,Inputfile, DBtitle, dbType))

def BLASTN_v29(outfile, infile, dbfile):
    import os
    import sys
    print "BLASTn (Sequence Read Match)"
    BLASTN_CMDLINE = "%sblastn -task blastn-short -out %s -query %s -db %s -outfmt 5 -num_alignments 250000"
    print BLASTN_CMDLINE, "^^^^^^^^^^^^^^^^^^^^^^^^^"
    status = os.system(BLASTN_CMDLINE % (UtilDir, outfile, infile, dbfile))
### NEED TO FINISH!!!


def BLASTN_v29template(outfile, infile, dbfile):
    import os
    import sys
    print "BLASTn (TEMPLATE SEARCH)"
    BLASTN_CMDLINE = "%sblastn -task blastn-short -out %s -query %s -db %s -outfmt 5 -num_alignments 1"
    print BLASTN_CMDLINE, "^^^^^^^^^^^^^^^^^^^^^^^^^"
    status = os.system(BLASTN_CMDLINE % (UtilDir, outfile, infile, dbfile))


def Analyze_Gal(ResultTableName, BLASTSeroDB):
    import MySQLdb
    import sys
    import os
    
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    CycleCounter = 0
    
    print "Searching for Template..."
    
    executelinePL = ("""SELECT * from FMD_%s""") % (ResultTableName)
    cursor.execute(executelinePL)

    Probe_List = cursor.fetchall()

    for Probe_Sel in Probe_List:

        ID_num = Probe_Sel[0]
        
        ID_name = Probe_Sel[1]
        Mean532 = Probe_Sel[2]
        Probe_Seq = Probe_Sel[3]
        #print Probe_Seq
        
        if int(float(Mean532)) >= int(Template_Threshold):  # Value error due to ".0" addition to the Mean532 value during conversion to integer. Fixed with int(float()).
            CycleCounter = CycleCounter + 1
            if CycleCounter == 5:                
                sys.stdout.write(".")
                CycleCounter = 0
                #try:
            Serotype_PrelimBLAST(ID_num, Probe_Seq, ResultTableName, BLASTSeroDB)
                #except:
                #print "error 171"
#continue

def Save_fasta(filename, title, sequence):
    fasta_file = open(filename, "w")
   # print "open 176 +++++"
    fasta_file.write(">" + title + "\n")
    for i in range(0, len(sequence), 72):
        fasta_file.write(sequence[i:i+72])
        fasta_file.write("\n")
    fasta_file.close()
    #print"close 176 -----"

### This BLAST step seeks out the best template ###

def BlastbyWindowDB(SeroDBinFile, Range, ResultTableName):
    from Bio import SeqIO
   
    SeroTempDB = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/testdb/TempTestDB.fasta"
    outputfileCON = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplatePREcon.fasta"

    openCON = open(outputfileCON, "w")    
    openCON.close()
    
    openSourcePRE = open(SeroDBinFile, "r") 
    print "open 197" 

    RangeGAP = Range*.8
    print Range, "<<<<<<<<<<< Range"
    RangeCount = int(4000/(RangeGAP))  #Changed from 10000
    RangeStart = 0
    RangeEnd = Range
        #try:
    for i in range(0,RangeCount):
        #try: 
            RangeStart = RangeStart + RangeGAP
            RangeEnd = RangeEnd + RangeGAP
            print RangeStart, RangeEnd
            
            if RangeStart > 10000:
                raise Exception("noobar")
            
            BuildWindowDB(SeroDBinFile, int(RangeStart), int(RangeEnd))
            
            
            Analyze_Gal(ResultTableName, "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_SelectedTemplateDB/FMDTemplateDB")
            #try:
                
            Get_Top_Scoring_Template_STEP1(ResultTableName)
            Get_Top_Scoring_Template_STEP2("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_SelectedTemplateDB/FMDTemplateDB", ResultTableName, RangeStart, RangeEnd, outputfileCON)
    
            
    openSourcePRE.close()        
    print "close 197"
    
def BuildWindowDB(SeroDBinFile,RangeStart, RangeEnd):
    from Bio import SeqIO
    if 1==1:
    #try
        SeroTempDB = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/testdb/TempTestDB.fasta"
          
        BlastWindow = open(SeroTempDB, "w")
        print "open 239"
        openSource = open(SeroDBinFile, "r")
        print "open241"
         
        for f in SeqIO.parse(openSource, "fasta"):
            PreOligo = f.seq.tostring()
            Oligo = PreOligo[RangeStart:RangeEnd]
            OligTitle = f.id[0:200]
            print OligTitle
            print Oligo
            
            OligTest = Oligo.count("-")
            if OligTest < (len(Oligo)*.85):
            
                BlastWindow.write(">" + OligTitle + "\n")
                for i in range(0, len(Oligo), 72):
                    BlastWindow.write(Oligo[i:i+72])
                    BlastWindow.write("\n")
     
       
        BlastWindow.close()
        print "close 239"
        openSource.close()
        print "close 241"   
          
        Format_SelectedTemplate_DB(SeroTempDB, "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_SelectedTemplateDB/FMDTemplateDB")


def Serotype_PrelimBLAST(ProbeID, ProbeSEQ, ResultTableName, Sero_DB):
     #print "RUNNING PRELIM BLAST..."
     from Bio import SeqIO
     from Bio.Blast import NCBIStandalone
     import MySQLdb

     PROBEID = str(ProbeID)
     PROBESEQ = ProbeSEQ
     PROBESEQlen = str(len(ProbeSEQ))
     AlignCNT = 0
     MaxAlignCNT = 1000
     COUNT = 0

     Save_fasta("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/ProbeBlastSeqPRELIM.fasta", PROBEID, PROBESEQ)
     
     BLASTN_v29template("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/CurrentFMDBlastPRELIM.xml","/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/ProbeBlastSeqPRELIM.fasta", Sero_DB)


     conn = MySQLdb.connect(host = HOSTlocal,
                                 user = USER,
                                 passwd = PASS,
                                 db = DB)  
     cursor = conn.cursor()

     result_handle = open("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/CurrentFMDBlastPRELIM.xml","r")
     
     from Bio.Blast import NCBIXML
     blast_records = NCBIXML.parse(result_handle)
     for blast_record in blast_records:
         for alignment in blast_record.alignments:
             for hsp in alignment.hsps:
                 COUNT = COUNT + 1
                 if AlignCNT < MaxAlignCNT:
                      AlignCNT = AlignCNT+1

                      SStart = (hsp.query_start)
                      SEnd = (hsp.query_end)
                     
                      PreGI = alignment.title
                      GIstart = 3                      
                      GIend = PreGI.find("|")
                      GIb = str(PreGI[GIstart:-1])
                     
                      GIbx = GIb.find("|")
                     
                      GI = str(GIb[0:GIbx])     
                      
                      COUNTs = str(COUNT)
                      IDENTITIES = str(hsp.identities)                      
         
                      Arguements = " '"+ COUNTs + "','" + GI + "','" + PROBEID + "','" + PROBESEQlen + "','" + IDENTITIES + "' "
                      EnterLine = "Prelim_Scoring_" + ResultTableName + " (Count, GI, Probe_ID, Probe_Length, Identities) VALUES (" + Arguements + ")"
                      ActLine = "INSERT INTO " + EnterLine
                     

                      cursor.execute(ActLine)
     
     cursor.close()
     conn.commit()
     conn.close()
     
def RunCAP(INFILE,OUTFILE):
    print "RUNNING CAP..."
    import os
    import MySQLdb
    from Bio import SeqIO
    from decimal import Decimal 
      
    #FASTACMD_CMDLINE = "c:\\CAP3\\CAP3 " + INFILE + " > " + OUTFILE + " -o 20 -p 80"  # CAP3 is too stringent, and it won't assemble all of the contigs for test cases.
    FASTACMD_CMDLINE = "wine /users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Tools/CAP.exe " + INFILE + " " + OUTFILE + " 20 80"
    print FASTACMD_CMDLINE   
    status = os.system(FASTACMD_CMDLINE)
    endStatus = os.system("/n") #attempting to send <CR> to end CAP
        
def Get_Top_Scoring_Template_STEP1(ResultTableName):    
    import MySQLdb
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    print "Scoring Templates..."

    executelineTS = ("""SELECT prelim_scoring_%s.GI, prelim_scoring_%s.Probe_ID,
                        Max(prelim_scoring_%s.Identities/prelim_scoring_%s.Probe_Length)
                        AS Ratio
                        FROM prelim_scoring_%s
                        GROUP BY prelim_scoring_%s.GI, prelim_scoring_%s.Probe_ID;
                        """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName)
    
    cursor.execute(executelineTS)
    TemplateScoring = cursor.fetchall()
    for TemplateScore in TemplateScoring:

        GI = str(TemplateScore[0])
        Probe_ID = str(TemplateScore[1])
        RATIOpre = str(TemplateScore[2])
        RATIO = RATIOpre[0:9]
        print RATIO 
        
        Arguements = " '" + GI + "','" + Probe_ID + "','" + RATIO + "' "
        EnterLine = "Scored_Probes_" + ResultTableName + " (GI, Probe_ID, Ratio) VALUES (" + Arguements + ")"
        ActLine = "INSERT INTO " + EnterLine
                
        cursor.execute(ActLine) 
   
    cursor.close()
    conn.commit()
    conn.close()

def Get_Top_Scoring_Template_STEP2(Sero_DB, ResultTableName, oStart,oEnd, filenameCON): 
    import MySQLdb
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    executelineTS2 = ("""SELECT Scored_Probes_%s.GI, Sum(Scored_Probes_%s.Ratio)
                        AS SumOfRatio
                        FROM Scored_Probes_%s
                        GROUP BY Scored_Probes_%s.GI
                        ORDER BY Sum(Scored_Probes_%s.Ratio) DESC
                        """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName)

    cursor.execute(executelineTS2)
    
    cursor.execute(executelineTS2)
    
    TemplateScoring2 = cursor.fetchone()
        #try:
    print TemplateScoring2
    Selected_Template = TemplateScoring2[0]
    
    Select_Template(Sero_DB, Selected_Template,oStart,oEnd, filenameCON)
        #except:
        #print "error 414"
        #print 'no template found...'
    cursor.close()
    conn.commit()
    conn.close()    

### Extracts Selected Template from Database, and creates a new BLAST database to build the consensus sequence ###

def Select_Template(Sero_DB,qGI, oStart, oEnd, filenameCON):
    import os
    import sys
    from Bio import SeqIO
    import MySQLdb

    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)    
    cursor = conn.cursor()    
   
    outputfile = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplate.fasta"
    #Deprecated 12/3/2014# FASTACMD_CMDLINE = "/Applications/bin/fastacmd -d %s -s %s -o %s"
    #status = os.system(FASTACMD_CMDLINE % (Sero_DB, qGI, outputfile))
    FASTACMDfile(qGI,Sero_DB,outputfile)
    
    if 1==1:
    #with open(outputfile) as handle:
  
        handle = open(outputfile)
        print "open 449"
        
        for seq_record in SeqIO.parse(handle, "fasta") :
            print " "
            print "*****************************  Selected Template  ******************************"
            print seq_record.description
            print "********************************************************************************"
            print " "
          
            oligSeq = seq_record.seq.tostring()
            print oStart
            print oEnd
            print oligSeq
           
            print 'SELOligoSequence >>>>>>>>>>>', oligSeq
            
            if oligSeq=="":
                print "Error 456 raised"
                raise Exception('foobar')
            
            Arguements = " '"+ str(oStart) + "','" + str(oEnd) + "','" + seq_record.description + "','" + oligSeq + "' "
            EnterLine = "WindowTemplate(StartPos, EndPos, Title, Sequence) VALUES (" + Arguements + ")"
            ActLine = "INSERT INTO " + EnterLine
            #print ActLine
            cursor.execute(ActLine)

            Title = seq_record.description
            
            #with open(filenameCON, "a") as _fileCONS:     
            fasta_fileCONS = open(filenameCON, "a")
            print "open 478"
            fasta_fileCONS.write(">" + Title + "\n")
            for i in range(0, len(oligSeq), 72):
                fasta_fileCONS.write(oligSeq[i:i+72])
                fasta_fileCONS.write("\n")
            fasta_fileCONS.close()
            print "close 478"

            
        Format_SelectedTemplate_DB(filenameCON,"/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_SelectedTemplateDB/FMDTemplateDB")  #<<< Switched from 'outputfile' to reflect window shifting
    
    
        handle.close()
        print "close 449"
        cursor.close()
        conn.commit()
        conn.close()
   
  

def Format_SelectedTemplate_DB(Inputfile, Outputfile):
    import os
    import sys
    print "RUNNING FORMAT_SELECTEDTEMPLATE_DB..."
    print "Creating Template Database..."
    
    #FORMATDB_CMDLINE = "formatdb -i %s -p F -o T -n c:\\FMDserotypingARRAY\\FMD_Selected_Template\\FMD_SelectedTemplateDB\\FMDTemplateDB"
    #FORMATDB_CMDLINE = "/Applications/bin/formatdb -i %s -p F -o T -n %s"
#status = os.system(FORMATDB_CMDLINE % (Inputfile, Outputfile))
    FormatDB(Inputfile, "F", Outputfile)
#MakeDB(Inputfile, "nucl", Outputfile)

############################################################## START CONSENSUS BUILD  #############################################################################

def Probe_Sort(ResultTableName, ProbeMFIThreshold):
    import MySQLdb

    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    dropline1 =("DROP TABLE IF EXISTS RESULTS_%s") % (ResultTableName)
    cursor.execute(dropline1)  
    executeline1 = ("""
            CREATE TABLE RESULTS_%s
            (
            Position VARCHAR(5),
            Probe_ID VARCHAR(5),
            Nuc_Score VARCHAR(10),
            Nucleotide VARCHAR(4))
            """) % (ResultTableName)
 
    cursor.execute (executeline1) 


    DeleteLine1 = """DELETE FROM RESULTS_SCORED_%s""" % (ResultTableName)
    DeleteLine2 = """DELETE FROM RESULTS_%s""" % (ResultTableName)
    DeleteLine3 = """DELETE FROM RESULTS_prescore_%s""" % (ResultTableName)
    DeleteLine4 = """DELETE FROM RESULTS_prescore2_%s""" % (ResultTableName)

    #cursor.execute(DeleteLine1)<< Stopped due to inclusion of MFI for final consensus
    #cursor.execute(DeleteLine2)
    cursor.execute(DeleteLine3)
    cursor.execute(DeleteLine4)
     
    #cursor.execute (executeline6)
    
    print "Scoring probes..."

    executelinePL = ("""SELECT * from FMD_%s""") % (ResultTableName)
    cursor.execute(executelinePL)
    print 'fetching probes'
    Probe_List = cursor.fetchall()

    for Probe_Sel in Probe_List:
        
        ID_num = Probe_Sel[0]
        ID_name = Probe_Sel[1]
        Mean532 = Probe_Sel[2]
        Probe_Seq = Probe_Sel[3]
        #print 'ID:::', ID_num, ' / ', Mean532
        
        if int(float(Mean532)) >= int(ProbeMFIThreshold):
            Serotype_BLAST(ID_num, Probe_Seq, Mean532, ResultTableName, ProbeMFIThreshold)
    print "DONE FETCHING PROBES..."
    cursor.close()
    conn.commit()
    conn.close()
       
def Serotype_BLAST(ProbeID, ProbeSEQ, ProbeScore, ResultTableName, ProbeMFIThreshold):
     from Bio import SeqIO
     from Bio.Blast import NCBIStandalone
     from Bio.Blast import NCBIXML
     import MySQLdb

     PROBEID = str(int(float(ProbeID)))
     PROBESEQ = ProbeSEQ     
     MaxAlignCNT = 1000
     AlignCNT = 0

     conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
     cursor = conn.cursor()
        
     Save_fasta("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/ProbeBlastSeq.fasta", PROBEID, PROBESEQ)

     Template_DB = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_FinalConsensusDB/FMDFinalConsensusDB"
     
     BLASTN_v29template("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/CurrentFMDBlast.xml","/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/ProbeBlastSeq.fasta", Template_DB)

     #print "BLASTING"
     result_handle = open("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/CurrentFMDBlast.xml","r")

     blast_records = NCBIXML.parse(result_handle)
     for blast_record in blast_records:
         for alignment in blast_record.alignments:
             for hsp in alignment.hsps:
                 if AlignCNT < MaxAlignCNT:
                      AlignCNT = AlignCNT+1

                      SStart = (hsp.query_start)
                      SEnd = (hsp.query_end)
                      
                      #print alignment.title
                      #print "Identities: " + str(hsp.identities)
                      #print "que: " + hsp.query +  str(hsp.query_start) + "::" + str(hsp.query_end)
                      #print "mat: " + hsp.match
                      #print "sub: " + hsp.sbjct + str(hsp.sbjct_start) + "::" + str(hsp.sbjct_end)
                      #print "Query Start NT:" + str(hsp.query_start)   
                      #print "Query Start NT:" + str(hsp.query_end)
                      #print "Subject Start NT:" + str(hsp.sbjct_start)
                      #print "Subject End NT:" + str(hsp.sbjct_end)
                      #print ""

                      preNUC_SCORE = str(float(ProbeScore)/len(hsp.sbjct))
                      NUC_SCORE = preNUC_SCORE[0:10]
                      
                      SubjectLength = len(hsp.query)
                      for NT in range(0,SubjectLength):
                          NT_ATCG = str(hsp.query[NT:NT+1])  
                          NT_Pos = str(hsp.sbjct_start+NT)

                          SERO = alignment.title
                                                   
                          Arguements = " '"+ NT_Pos + "','" + PROBEID + "','" + NUC_SCORE + "','" + NT_ATCG + "' "
                          EnterLine = "RESULTS_" + ResultTableName  + " (Position, Probe_ID, Nuc_Score, Nucleotide) VALUES (" + Arguements + ")"
                          ActLine = "INSERT INTO " + EnterLine
                          
                          cursor.execute(ActLine)
                          #print 'done'
     
     cursor.close()
     conn.commit()
     conn.close()
    # print 'finished'
def Score_Consensus(ResultTableName, ProbeMFIThreshold):
    import MySQLdb
    print "Score Consensus..."
    conn1 = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor1 = conn1.cursor()

    print "Building Consensus..."


    

    executelineSC = ("""SELECT results_%s.Position,
                        Sum(results_%s.Nuc_Score)
                        AS SumOfNuc_Score, results_%s.Nucleotide
                        FROM results_%s
                        GROUP BY results_%s.Position, results_%s.Nucleotide;
                        """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName)
    
    cursor1.execute(executelineSC)
    Scoring_List = cursor1.fetchall()
    for Probe_Scores in Scoring_List:

        NT_Pos = str(Probe_Scores[0])
        preSUMOFSCORE = str(Probe_Scores[1])
        SUMOFSCORE = preSUMOFSCORE[0:10]
        NT_ATCG = Probe_Scores[2]

        #print NT_Pos
        #print SUMOFSCORE
        #print NT_ATCG

        Arguements = " '"+ NT_Pos + "','" + SUMOFSCORE + "','" + NT_ATCG + "' "
        EnterLine = "RESULTS_preScore_" + ResultTableName + " (Position, SumOfScore, Nucleotide) VALUES (" + Arguements + ")"
        ActLine = "INSERT INTO " + EnterLine
        cursor1.execute(ActLine)

    cursor1.close()
    conn1.commit()
    conn1.close()

    conn2 = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor2 = conn2.cursor()

    executelineSC2 = ("""SELECT RESULTS_preScore_%s.Position,
                        Max(RESULTS_preScore_%s.SumOfScore)
                        AS MaxOfSumOfNuc_Score
                        FROM RESULTS_preScore_%s
                        GROUP BY RESULTS_preScore_%s.Position;
                        """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName)
    
    cursor2.execute(executelineSC2)
    Scoring_List2 = cursor2.fetchall()
    for Probe_Scores2 in Scoring_List2:

        NT_PosMAX = str(Probe_Scores2[0])
        MAXOFSCORE = str(Probe_Scores2[1])

        #print NT_PosMAX
        #print MAXOFSCORE
        
        Arguements2 = " '"+ NT_PosMAX + "','" + MAXOFSCORE +  "' "
        EnterLine2 = "RESULTS_preScore2_" + ResultTableName + " (Position, MaxOfScore) VALUES (" + Arguements2 + ")"
        ActLine2 = "INSERT INTO " + EnterLine2
        cursor2.execute(ActLine2)

    cursor2.close()
    conn2.commit()
    conn2.close()

    conn3 = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor3 = conn3.cursor()

    executelineSC3 = ("""SELECT RESULTS_preScore_%s.Position, RESULTS_preScore_%s.Nucleotide
                        FROM RESULTS_preScore2_%s
                        INNER JOIN RESULTS_preScore_%s
                        ON (RESULTS_preScore2_%s.MaxOfScore = RESULTS_preScore_%s.SumOfScore)
                        AND (RESULTS_preScore2_%s.Position = RESULTS_preScore_%s.Position);
                        """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName)
    
    cursor3.execute(executelineSC3)
    Scoring_List3 = cursor3.fetchall()
    for Probe_Scores3 in Scoring_List3:

        NT_PosFINAL = str(Probe_Scores3[0])
        NucleotideFINAL = str(Probe_Scores3[1])

        print NT_PosFINAL
        print NucleotideFINAL

        Arguements3 = " '"+ NT_PosFINAL + "','" + NucleotideFINAL + "','" + str(ProbeMFIThreshold) + "' " #<< Added MFI for building final Consensus
        EnterLine3 = "RESULTS_Scored_" + ResultTableName + " (Position, Nucleotide, MFI) VALUES (" + Arguements3 + ")" #<< Added MFI for building final Consensus
        ActLine3 = "INSERT INTO " + EnterLine3
        cursor3.execute(ActLine3)

    cursor3.close()
    conn3.commit()
    conn3.close()
    

def Build_Consensus(ResultTableName, OutputFileName, ResultTitle, ProbeMFIThreshold):
    import MySQLdb
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()


    executelineCON = ("""SELECT * from RESULTS_Scored_%s""") % ResultTableName
    cursor.execute(executelineCON)

    FileOUT_Name = """/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/%s.FASTA""" % (OutputFileName)
    
    Result_File = open(FileOUT_Name, "a")
    print "open 778"
    Result_File.write(">" + ResultTitle + "_" + str(ProbeMFIThreshold) + "_CONSENSUS"  + "\n")
    
    for i in range (1, 4000):         
        ConsensusSelect = ("""SELECT RESULTS_Scored_%s.Position, RESULTS_Scored_%s.Nucleotide
                            FROM RESULTS_Scored_%s
                            GROUP BY RESULTS_Scored_%s.Position, RESULTS_Scored_%s.Nucleotide, RESULTS_Scored_%s.MFI
                            HAVING (((RESULTS_Scored_%s.Position)=%s) AND ((RESULTS_Scored_%s.MFI)=%s));""") % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, i, ResultTableName, ProbeMFIThreshold)
        
        cursor.execute(ConsensusSelect)        
        Consensus_Hit = cursor.fetchone()

        if Consensus_Hit == None:
            NT_Sel = "N"
        else:
            NT_Sel = Consensus_Hit[1]

        print str(i) + ".  " + NT_Sel
        Result_File.write(NT_Sel)
    Result_File.write("\n")
    Result_File.close()
    print "close 778"


    
def Parse_GPR_file(GPR_filename):
    import MySQLdb
    
    conn = MySQLdb.connect(host = HOSTremote,
                       user = USER,
                       passwd = PASS,
                       db = DB)
    cursor = conn.cursor()

#cursor.execute("DELETE FROM FMD_GPR")
#cursor.execute("DELETE FROM FMD_GPRLONG")
#cursor.execute("DELETE FROM FMD_MIDb")
#cursor.execute("DELETE FROM FMD_GPRXLONG")

    print "Parsing GPR file..." 

    GPR_file = open(GPR_filename, "r")
    Clean_GPR_file = open("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/CleanGPR.txt", "w")
    LineCounter = 0
    STARTGPR = 'FALSE'
    countertest =0
    for line in GPR_file:
        if countertest <=35:
            print line[0:5]
            countertest = countertest+1
        #STARTGPR = STARTGPR + 1
        if line[1:6] == "Block":
            print "YEAH!"
            LineCounter = LineCounter + 1
            STARTGPR = 'TRUE' 
       
            continue
        
        if LineCounter > 0 and STARTGPR == 'TRUE':
            #print line
            LenLine = len(line)-1
            Clean_GPR_file.write(line[0:LenLine])
            Clean_GPR_file.write("\t")
            Clean_GPR_file.write("\n")
    Clean_GPR_file.close()
    GPR_file.close()
  
    EnterLine = """LOAD DATA LOCAL INFILE '/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/CleanGPR.txt'
                  INTO TABLE FMD_GPR"""
    
    cursor.execute(EnterLine)                    

    cursor.close()
    conn.commit()
    conn.close()

def Build_Analysis_File(ResultTableName):
    import MySQLdb
    
    conn = MySQLdb.connect(host = HOSTremote,
                       user = USER,
                       passwd = PASS,
                       db = DB)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM FMD_GPR")
    GALFILE = cursor.fetchall()

    for GFline in GALFILE:
        
        ID_NUM = str(GFline[40])  #Modified this from 40 to 30.  Possible difference between GAL and GPR files??? 9/13/13
        ID_NAME = str(GFline[3])
        MEAN532 = str(GFline[9])
        PROBE_SEQ = str(GFline[44]) #Modified this from 44 to 34.  Possible difference between GAL and GPR files??? 9/13/13
        
        if len(PROBE_SEQ) > 2:
            
            ArguementsGF = " '"+ ID_NUM + "','" + ID_NAME + "','" + MEAN532 + "','" + PROBE_SEQ +  "' "
            EnterLineGF = "FMD_" + ResultTableName + " (ID_Num, ID_Name, Mean532, Probe_Seq) VALUES (" + ArguementsGF + ")"
            ActLineGF = "INSERT INTO " + EnterLineGF
            cursor.execute(ActLineGF)   

    cursor.close()
    conn.commit()
    conn.close() 

def FinalConsensusBuild(ResultTableName, Threshold):
    
    Probe_Sort(ResultTableName,Threshold)
    Score_Consensus(ResultTableName, Threshold)
    Build_Consensus(ResultTableName, ResultTableName, ResultTableName, Threshold)

def Align_Results(OutputFileName):
    import os
    
    FileIN_Name = """/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/%s.FASTA""" % (OutputFileName)
    FileOUT_ALN = """/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/%s.ALN""" % (OutputFileName)
    print FileIN_Name
    print FileOUT_ALN
    
    from Bio.Clustalw import MultipleAlignCL
    from Bio import Clustalw

    cline = MultipleAlignCL(os.path.join(os.curdir, FileIN_Name))
    cline.set_output(FileOUT_ALN)
    
    alignment = Clustalw.do_alignment(cline)

    cline.close()
    
def Build_FINAL_Consensus(ResultTableName):
    import MySQLdb
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()

    #dropline7 =("DROP TABLE IF EXISTS RESULTS_PREconsensus_%s")%(ResultTableName)  #<<< Temporary Addition for scoring
    #cursor.execute(dropline7)
    #executeline7 = ("""
    #        CREATE TABLE RESULTS_PREconsensus_%s
    #        (
    #        Position VARCHAR(5),
    #        MaxMFI VARCHAR(10))            
    #        """)% (ResultTableName) #<< Added table for building final Consensus

   # cursor.execute (executeline7)
    
    PreConsensus = ("""SELECT results_scored_%s.Position,
                        Max(results_scored_%s.MFI) AS MaxOfMFI
                        FROM results_scored_%s
                        GROUP BY results_scored_%s.Position; """)% (ResultTableName, ResultTableName, ResultTableName, ResultTableName)

    cursor.execute(PreConsensus)
    PreCONS = cursor.fetchall()
    for precon in PreCONS:
        PositionPC = precon[0]
        MaxMFIPC = precon[1]
        
        ArguementsPC = " '"+ str(PositionPC) + "','" + str(MaxMFIPC) +  "' "
        EnterLinePC = "RESULTS_PREconsensus_" + ResultTableName + " (Position, MaxMFI) VALUES (" + ArguementsPC + ")"
        ActLinePC = "INSERT INTO " + EnterLinePC
        cursor.execute(ActLinePC)   

    cursor.close()
    conn.commit()
    conn.close()
    
    Build_FINAL_Consensus2(ResultTableName)

def Build_FINAL_Consensus2(ResultTableName):
    import MySQLdb
    conn = MySQLdb.connect(host = HOSTlocal,
                           user = USER,
                           passwd = PASS,
                           db = DB)
    cursor = conn.cursor()



    FileOUT_Name = """/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/Consensus_Results/%s.FASTA""" % (ResultTableName)
    
    Result_File = open(FileOUT_Name, "a")
    Result_File.write(">" + ResultTableName + "_" + "Final_CONSENSUS"  + "\n")
    
    for i in range (1, 4000):         
        ConsensusSelect = ("""SELECT results_scored_%s.Position, results_scored_%s.Nucleotide
                                FROM RESULTS_PREconsensus_%s
                                INNER JOIN results_scored_%s
                                ON (RESULTS_PREconsensus_%s.MaxMFI = results_scored_%s.MFI)
                                AND (RESULTS_PREconsensus_%s.Position = results_scored_%s.Position)
                                WHERE (((results_scored_%s.Position)=%s));
                                """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, i)


        cursor.execute(ConsensusSelect)        
        Consensus_Hit = cursor.fetchone()
        #print Consensus_Hit[0]
        if Consensus_Hit == None:
            NT_Sel = "N"
        else:
            NT_Sel = Consensus_Hit[1]

        print str(i) + ".  " + NT_Sel
        Result_File.write(NT_Sel)
    Result_File.write("\n>" + ResultTableName + "_" + "Final_CONSENSUS_SCORE"  + "\n")    
    for j in range (1, 4000):         
        ConsensusSelectB = ("""SELECT results_scored_%s.Position, results_scored_%s.Nucleotide, results_scored_%s.MFI
                                FROM RESULTS_PREconsensus_%s
                                INNER JOIN results_scored_%s
                                ON (RESULTS_PREconsensus_%s.MaxMFI = results_scored_%s.MFI)
                                AND (RESULTS_PREconsensus_%s.Position = results_scored_%s.Position)
                                WHERE (((results_scored_%s.Position)=%s));
                                """) % (ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, ResultTableName, j)

        
        cursor.execute(ConsensusSelectB)        
        Consensus_Score = cursor.fetchone()
        NT_Score = Consensus_Score[1]
        if NT_Score == "1562":
            NT_Code = 'F'
        else:
            if NT_Score == "3125":
                NT_Code = 'E'
            else:
                if NT_Score == "6250":
                    NT_Code = 'D'
                else:
                    if NT_Score == "12500":
                        NT_Code = 'C'
                    else:
                        if NT_Score == "25000":
                            NT_Code = 'B'
                        else:
                                                
                            if NT_Score == "50000":
                                NT_Code = 'A'
                            else:
                                NT_Code = 'X'
        print str(j) + ".  " + NT_Code
        Result_File.write(NT_Code)
    
    Result_File.write("\n")
    Result_File.close()
    


#####################################################
#           Main Program
#####################################################

# --DATABASE VARIABLE--
HOSTlocal = "localhost"
HOSTremote = "localhost"
DB = "fmdserotyping"
USER = "root"
PASS = "Piobmare,774"
ComputerID = "rbarrette"
AlignCNT = 0
WorkingDir = "/users/rwbarrettemac/bioinformatics/pythonfolders/TorrentFiles/"
UtilDir = "/usr/local/ncbi/blast/bin/"


#FMD serotyping Reference Genomes for blast database#
#FMDserotypeDB = "C:\\FMDserotypingARRAY\\FMD_NewRef_DB\\FMD_CompleteGenomeDB"
FMDserotypeDB = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/testdb/FMDrefTest"

#Probe Intensity Cutoff
Template_Threshold = 10000


def RunFMDAnalysis(folder, fileR):#

    worklist = folder + fileR
    worklistfile = open(worklist,"r")
        #try:
    for line in worklistfile:
            print "---> Starting analysis on " + line#

        
            ResultTableNamePre = line.split()
            
            ResultTableName = ResultTableNamePre[0]#
            print ResultTableNamePre
            print ResultTableName
            
            Create_Results_Table(ResultTableName, ResultTableName, "TRUE") #Added 9/13/13 to reinitialize tables
        
            gprfile = str(folder + ResultTableName + ".gpr")
            print gprfile
            Parse_GPR_file(gprfile)    # Added 9/13/2013 to deal with script not handling GPR files
            Build_Analysis_File(ResultTableName) # Added 9/13/2013 to deal with script not handling GPR files to create intermidate FMD_xxx
        
        #BlastbyWindowDB("C:\\FMDserotypingARRAY\\FMDcompleteGenomeALN.fasta", 500, ResultTableName)
            BlastbyWindowDB("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDdata/fmdcgallALN_NEWtrunc.fasta", 500, ResultTableName)
            ContigIN = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingArray/FMD_Selected_Template/Selected_FMDTemplatePREcon.fasta.cap.contigs"
        ##ContigOUTold = "c:\\FMDserotypingArray\\consensus_results\\" + ResultTableName + "ContigOUT.fasta"
            ContigOUT = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingArray/FMD_Selected_Template/Selected_FMDTemplateContigOUT.fasta"
            ConsensusOUT = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingArray/FMD_Selected_Template/FMD_FinalConsensusDB/FMD_FinalConsensus.fasta"
        
            RunCAP("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplatePREcon.fasta","/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplatePOSTconX.fasta")
            # ADDED BELOW
            #  ParseContig(ContigOUT, ConsensusOUT)
            #Format_SelectedTemplate_DB(ConsensusOUT, "/users/rogerbarrette/documents/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_FinalConsensusDB/FMDFinalConsensusDB")

#NewFMDconsBuild(ResultTableName)

#           Build_FINAL_Consensus(ResultTableName)
###MOVE### RunCAP3("c:\\FMDserotypingARRAY\\FMD_Selected_Template\\Selected_FMDTemplatePREcon.fasta","c:\\FMDserotypingARRAY\\FMD_Selected_Template\\Selected_FMDTemplatePOSTcon")
        
        ###################  << Too many files open at this point....not sure what is still alive
        ###MOVE### CompileContigs('1000', ContigIN, ContigOUT)

        ###MOVE### Format_SelectedTemplate_DB(ContigOUT)
        
        ###MOVE### NewFMDconsBuild(ResultTableName)
        
     #   ###MOVE### Build_FINAL_Consensus(ResultTableName)
     ###except:
     ###  print "Done..."
    worklistfile.close()
    
def RunFMDAnalysisCONSENS(folder, fileR):

    worklist = folder + fileR
    worklistfile = open(worklist,"r")
    for line in worklistfile:
        #try:
            print "---> Starting Final Consensus Assembly " + line

            
            ResultTableNamePre = line.split()
            ResultTableName = ResultTableNamePre[0]


           
            #ContigIN = "c:\\FMDserotypingArray\\FMD_Selected_Template\\Selected_FMDTemplatePREcon.fasta.cap.contigs"  Deprecated?
            #ContigOUT = "c:\\FMDserotypingArray\\FMD_Selected_Template\\Selected_FMDTemplateContigOUT.fasta" # Changed due to issues with CAP3 inability to produce the consensus
            ContigOUT = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingArray/FMD_Selected_Template/Selected_FMDTemplatePOSTconX.fasta"
            ConsensusOUT = "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingArray/FMD_Selected_Template/FMD_FinalConsensusDB/FMD_FinalConsensus.fasta"
            
            
            #RunCAP("/users/rogerbarrette/documents/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplatePREcon.fasta","/users/rogerbarrette/documents/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/Selected_FMDTemplatePOSTconX.fasta")
            
            ###################  << Too many files open at this point....not sure what is still alive
            # DEPRECATED:::CompileContigs('1000', ContigIN, ContigOUT)
          
            # DEPRECATED:::Format_SelectedTemplate_DB(ContigOUT, "c:\\FMDserotypingARRAY\\FMD_Selected_Template\\FMD_FinalConsensusDB\\FMDFinalConsensusDB")
            
            ParseContig(ContigOUT, ConsensusOUT)
            Format_SelectedTemplate_DB(ConsensusOUT, "/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDserotypingARRAY/FMD_Selected_Template/FMD_FinalConsensusDB/FMDFinalConsensusDB")
            
            NewFMDconsBuild(ResultTableName)
            
            Build_FINAL_Consensus(ResultTableName)
                #except:
                #continue
    worklistfile.close()

def ParseContig(ContigIN, FastaOUT):
    import sys
    import MySQLdb
    from Bio import SeqIO 
    from decimal import Decimal
    
    fasta_file = open(FastaOUT, "w") #
    handle = open(ContigIN)

    for seq_record in SeqIO.parse(handle, "fasta"):
        title = seq_record.description 
        sequence = seq_record.seq.tostring()
        SeqLEN = len(seq_record.seq.tostring())
        
        if title == "Contig-0":
            print "TTTTTTTTTTTTTTRRRRRRRRRRRUUUUUUUUUUUUUUUEEEEEEEEEEEEEE!!!!"
            print title, sequence
            fasta_file.write(">" + "FMD_Template" + "\n")
            for i in range(0, len(sequence), 72):
                fasta_file.write(sequence[i:i+72])
                fasta_file.write("\n")
        if title == "Contig-1":
            raise "ERROR: too many contigs"
    fasta_file.close()
    handle.close()
   
    
        
def CompileContigs(qGI,ContigIN, ContigOUT):    
    import sys
    import MySQLdb
    from Bio import SeqIO 
    from decimal import Decimal
    
    fasta_file = open(ContigOUT, "w")
    print "open 1048"
    handle = open(ContigIN)
    print "open 1050"
    for seq_record in SeqIO.parse(handle, "fasta"):
        title = seq_record.description 
        sequence = seq_record.seq.tostring()
        SeqLEN = len(seq_record.seq.tostring())
        NEWtitle = qGI + "_" + title
        NEWtitleB = "SelectedTemplate"
        print title, sequence
            
        fasta_file.write(">" + NEWtitleB + "\n")
        for i in range(0, len(sequence), 72):
            fasta_file.write(sequence[i:i+72])
            fasta_file.write("<cr>")
    fasta_file.close()   
    print "close 1048"
    handle.close()
    print "close 1050"
    
def NewFMDconsBuild(ResultTableName):

       ##################  BUILD CONSENSUS  ###################

        FinalConsensusBuild(ResultTableName, 1562)
        FinalConsensusBuild(ResultTableName, 3125)
        FinalConsensusBuild(ResultTableName, 6250)
        FinalConsensusBuild(ResultTableName, 12500)
        FinalConsensusBuild(ResultTableName, 25000)
        FinalConsensusBuild(ResultTableName, 50000)
    #########################################################


RunFMDAnalysis("/users/rwbarrettemac/bioinformatics/pythonfolders/FMDanalysisScript/FMDdata/9Jun10/","worklist.txt")

        
        