#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define		NUM_BITS	4

#define CATCH_ERROR(error) {{if (error) { fprintf(stderr, "%s\n", processError(error)); exit(EXIT_FAILURE); }}}
#ifndef MIN
	#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif


typedef struct {
	uint size;
	uint numEntries; 
	uint *h_reference;
	uint *d_reference;
	unsigned char *char_reference;
} ref_t;

uint base2number(unsigned char base)
{
	switch(base)
	{
    	case 'A':
    	case 'a':
    	    return(0);
    	case 'C':
    	case 'c':
    	    return(1 << 28);
    	case 'G':
    	case 'g':
    	    return(2 << 28);
    	case 'T':
    	case 't':
    	    return(3 << 28);
    	default :
    	    return(4 << 28);
	}
}

//basesPerWord tiene que modificarse por si el tamaño de la referencia no es multiple del tamaño de la entrada
void char2bin(unsigned char *reference, uint *binaryRef, uint pos, uint basesPerWord)
{
	int i; 
	uint indexBase;
	uint bitmap	= 0x0;

	for(i = 0; i < basesPerWord; i++){
		indexBase = base2number(reference[pos + i]);
		bitmap >>= NUM_BITS;
		bitmap |= indexBase;
	}

	binaryRef[pos / basesPerWord] = bitmap;
}

int transformReference(void *reference)
{
	ref_t *ref = (ref_t *) reference;
	uint pos;
	uint basesPerWord = 32 / NUM_BITS;	

	printf("TR: size: %u, numentries: %u\n", ref->size, ref->numEntries);

	ref->h_reference = (uint*) malloc(ref->numEntries * sizeof(uint));
	for(pos = 0; pos < ref->size; pos += basesPerWord){
		char2bin(ref->char_reference, ref->h_reference, pos, basesPerWord);
	}
	
	return (0);	
}

int loadRef(const unsigned char *fn, void **reference)
{      
	ref_t *ref = (ref_t *) reference;
	FILE *fp = NULL;
	unsigned char cadena[256];
	uint numCharacters = 0;
	uint sizeFile = 0;
	uint cleidos = 0;
	uint pos = 0;
	uint basesPerWord = 32 / NUM_BITS;	

	int i = 0;

	fp = fopen(fn, "rb");
	if (fp==NULL) return (30);

	fseek(fp, 0L, SEEK_END);
	sizeFile = ftell(fp);
	rewind(fp);

	printf("sizeFile: %u\n", sizeFile);

	ref->char_reference = (unsigned char*) malloc(sizeFile * sizeof(char));
	if (ref==NULL) return (31);
	
	if ((fgets(cadena, 256, fp) == NULL) || (cadena[0] != '>')) 
		return (32);

	while((!feof(fp)) && (fgets(cadena, 256, fp) != NULL)){
		if (cadena[0] != '>'){
			cleidos = strlen(cadena);
			if(cleidos) cleidos--;
			memcpy((ref->char_reference + pos), cadena, cleidos);
			pos +=  cleidos;
		}
	}

	fclose(fp);

	ref->numEntries = (pos / basesPerWord) + ((pos % basesPerWord) ? 1 : 0);
	ref->size = pos;

	printf("size: %u, numentries: %u\n", ref->size, ref->numEntries);

	return (0);
}

int saveRef(const char *fn, void *reference)
{
    ref_t *ref = (ref_t *) reference;

    char fmiFileOut[512];
    FILE *fp=NULL;
    uint i, error;

    sprintf(fmiFileOut, "%s.%u.4bits.ref", fn, ref->size);
    
    fp = fopen(fmiFileOut, "wb");
    if (fp==NULL) return (8);

    fwrite(&ref->numEntries, sizeof(uint), 1, fp);
    fwrite(&ref->size, sizeof(uint), 1, fp);

    fwrite(ref->h_reference, sizeof(uint), ref->numEntries, fp);
    fclose(fp);
    return (0);
}

int freeReference(void *reference)
{   
    ref_t *ref = (ref_t *) reference;  

    if(ref->h_reference != NULL){
        free(ref->h_reference);
        ref->h_reference=NULL;
    }

    if(ref->char_reference != NULL){
        free(ref->char_reference);
        ref->char_reference=NULL;
    }   
	
    return(0);
}

char *processError(int e){ 
    switch(e) {
        case 0:  return "No error"; break; 
        case 30: return "Cannot open reference file"; break;
        case 31: return "Cannot allocate reference"; break;
        case 32: return "Reference file isn't multifasta format"; break;
        case 37: return "Cannot open reference file on write mode"; break;
        case 42: return "Cannot open queries file"; break;
        case 43: return "Cannot allocate queries"; break;
        case 45: return "Cannot allocate results"; break;
        case 47: return "Cannot open results file for save intervals"; break;
        case 48: return "Cannot open results file for load intervals"; break;
        case 99: return "Not implemented"; break;
        default: return "Unknown error";
    }   
}

int initReference(void **reference)
{
    ref_t *ref = (ref_t *) malloc(sizeof(ref_t));
	ref->size = 0;
	ref->numEntries = 0;

	ref->h_reference = NULL;
	ref->d_reference = NULL;
	ref->char_reference = NULL;

    (*reference) = ref;
    return (0);
}

int main(int argc, char *argv[])
{
    void *reference;
    unsigned char *refFile = argv[1];
    int error;

	error = initReference(&reference);    
	CATCH_ERROR(error);

	error = loadRef(refFile, reference);
    CATCH_ERROR(error);
	
	error = transformReference(reference);
    CATCH_ERROR(error);

    error = saveRef(refFile, reference);
    CATCH_ERROR(error);
    
    error = freeReference(reference);
    CATCH_ERROR(error);

    return (0);
}

