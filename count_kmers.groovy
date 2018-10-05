def cli = new CliBuilder(usage:'kmer_dig [options] message')

cli._(longOpt:'inputFileName', args:1, argName:'inputFileName', 'repertoire file, with columns referring to aaSeqCDR3, Vgene, cloneFraction')
cli.v(longOpt: 'Vtable', args:1, argName:'v','table with V genes, optional')
cli.d(longOpt: 'discont', args:1, argName: 'd','do we count discontinious k-mers, default = false')
cli.f(longOpt: 'useFreq',args:1, argName:'f' ,'do we use clone frequencies or clone numbers, default = false')
cli.l(longOpt: 'cutLeft',args:1, argName: 'l', 'number of amino acids to cut from the left side of CDR3, default = 3')
cli.r(longOpt: 'cutRight',args:1, argName: 'r' ,'number of amino acids to cut from the right side of CDR3,default = 3')
cli.lmax(longOpt: 'kmerLengthMax',args:1, argName: 'lmax','max length of the k-mer, default = 3')
cli.lmin(longOpt: 'kmerLengthMin',args:1, argName: 'lmin','min length of the k-mer, default = 3')

cli.h('print this message')

def options = cli.parse(args)

if (options.h){
	cli.usage()
	System.exit(0)
} else {

def inputFileName = options.arguments()[0]
def Vtable = options.v
def discont = options.d.toBoolean()
def useFreq = options.f.toBoolean()
def cutLeft = options.l ? options.l.toInteger() : 3    // ternary operator (if) to set default of cutLeft to3
def cutRight = options.r ? options.r.toInteger(): 3    // elvis operator (if) to set default of cutRight to 3
def lmax = options.lmax ? options.lmax.toInteger() : 4    // elvis operator (if) to set default of max k-mer length to 4
def lmin = options.lmin ? options.lmin.toInteger() : 3    // elvis operator (if) to set default of max k-mer length to 4

if (lmin<3) System.err << "\nwarning! min k-mer length is too short (less than 3). results for such short k-mers may be doubtful\n"

//println("vtable " +Vtable)
//println("disc "+ discont)
//println("freq " + useFreq)
//println("left " + cutLeft)
//println(cutRight)
//println(lmax)
//println(lmin)
//println("filename " + inputFileName)



def firstLine = true

//read V genes similarity table

def Vref = [:]
new File(Vtable).splitEachLine(":") {v ->
	Vref.put(v[0],v[1].split(",").sort())
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

	aaCDR3Col=seq.findIndexOf{it in ["aaSeqCDR3","CDR3aa", "CDR3.amino.acid.sequence"]}
	if(aaCDR3Col==-1) {
		System.err << "Could not find column with aa CDR3 sequences - no column named aaSeqCDR3 or CDR3aa\n"
		System.exit(0)
	}

	VCol=seq.findIndexOf{it in ["bestVGene","V", "V.gene"]}
        if(VCol==-1) {
                System.err << "Could not find column with V gene annotation - no column named bestVGene or V\n"
                System.exit(0)
        }

	if (useFreq) {
		freqCol=seq.findIndexOf{it in ["cloneFraction","frequency", "Read.proportion"]}
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

	
	if(seq[aaCDR3Col].length() > cutLeft+cutRight+1) {
		seq[aaCDR3Col]=seq[aaCDR3Col].substring(cutLeft,seq[aaCDR3Col].length()-cutRight);

		V = seq[VCol].split(",");
		for (k=lmin;k<=lmax;k++){
			currentKmers =[]
			for (i=0;i<seq[aaCDR3Col].length()-k+1;i++){
				kmer = seq[aaCDR3Col].substring(i,i+k);
				if (currentKmers.contains(kmer)) continue
				currentKmers << kmer
				Vmap = kmerMap.get(kmer)
				if (useFreq) {freq = Double.parseDouble(seq[freqCol])
					} else { freq = 1}
				V.each{v->
						if (Vtable != false){
							if (Vmap.containsKey((v))) {
                                                                kmerMap.get(kmer) << [(Vref.get(v).join('/')) : (Vmap.get((v)) + 1*freq)]
                                                        } else {
                                                                kmerMap.get(kmer)<< [(Vref.get(v).join('/')) :1*freq]
                                                        }

						} else {
								
							if (Vmap.containsKey((v))) {
								kmerMap.get(kmer) << [(v) : (Vmap.get((v)) + 1*freq)]
							} else {
								kmerMap.get(kmer)<< [(v):1*freq]
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
}
