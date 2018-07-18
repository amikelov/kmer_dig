
def inputFileName = args[0]
def Vtable = args[1]
def discont = args[2].toBoolean()
def useFreq = args[3].toBoolean()

assert args.length == 4

def firstLine = true

//read V genes similarity table

def Vref = [:]
new File(Vtable).splitEachLine(":") {v ->
	Vref.put(v[0],v[1].split(","))
}

def freq
def kmerMap=[:].withDefault{[:].withDefault{0}}
def kmer
int fileLength = 0
def V
def Vmap
def currentKmers = []

def aaCDR3Col = 0
def VCol = 0
def freqCol= 0

new File(inputFileName).splitEachLine("\t") {seq ->
   if(firstLine) {
	
	//find all needed columns

	aaCDR3Col=seq.findIndexOf{it in ["aaSeqCDR3","CDR3aa"]}
	if(aaCDR3Col==-1) {
		System.err << "Could not find column with aa CDR3 sequences - no column named aaSeqCDR3 or CDR3aa\n"
		System.exit(0)
	}

	VCol=seq.findIndexOf{it in ["bestVGene","V"]}
        if(VCol==-1) {
                System.err << "Could not find column with V gene annotation - no column named bestVGene or V\n"
                System.exit(0)
        }

	if (useFreq) {
		freqCol=seq.findIndexOf{it in ["cloneFraction","frequency"]}
		 if(freqCol==-1) {
			System.err << "Could not find column with V gene annotation - no column named cloneFraction or frequency\n"
			System.exit(0)
		}
	}

	firstLine= false
   } else {
	
	// UI counter, sends message each 10000 clones
	fileLength++
	if (!(fileLength % 100000)){
		System.err << fileLength.intdiv(100000)
		System.err << "\n"
	}

	
	if(seq[aaCDR3Col].length() > 7) {
		seq[aaCDR3Col]=seq[aaCDR3Col].substring(3,seq[aaCDR3Col].length()-3);

		V = seq[VCol].split(",");
		for (k=3;k<5;k++){
			currentKmers =[]
			for (i=0;i<seq[aaCDR3Col].length()-k+1;i++){
				kmer = seq[aaCDR3Col].substring(i,i+k);
				if (currentKmers.contains(kmer)) continue
				currentKmers << kmer
				Vmap = kmerMap.get(kmer)
				if (useFreq) {freq = Double.parseDouble(seq[freqCol])
					} else { freq = 1}
				V.each{v->
					Vref.get(v).each {entry ->		
						if (Vmap.containsKey((entry))) {
							kmerMap.get(kmer) << [(entry) : (Vmap.get((entry)) + 1*freq)]
						} else {
							kmerMap.get(kmer)<< [(entry):1*freq]
						}
					}
				}
			}
		}
	}
   }
}

System.err << "\nDone with continious, going on\n"

if (discont) {
def kmerDiscontMap = [:]
def kmerD
def regex
def nmatch=0
def result = [:].withDefault{0}
def maps =[]
def usedKeys = []
for (km in kmerMap.keySet()){
	for(l=0;l<km.length()-2;l++){
		kmerD = km.substring(0,l+1) + "." + km.substring(l+2,km.length())
		if (!usedKeys.contains(kmerD) && (kmerD.length()>3) ) {
			usedKeys << kmerD
			regex = ~"$kmerD"
			result = [:].withDefault{0}
			maps=[]
			kmerMap.keySet().grep(regex).each{key ->
				maps<<kmerMap.get(key)
			}		 
			maps.collectMany{it.entrySet()}.each{ result[it.key] += it.value }
			kmerDiscontMap<< ["${kmerD}":result]				
		}	
	}
}

kmerMap<<kmerDiscontMap
}

println(fileLength)
BufferedWriter log = new BufferedWriter(new OutputStreamWriter(System.out));
kmerMap.each{k,v ->
	v.each{Vgene,n ->
		log.write("${k}"+"\t"+"${Vgene}" +"\t"+n +"\n")
	
	}
}
log.flush()
