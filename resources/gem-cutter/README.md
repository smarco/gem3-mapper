# GEM-Cutter

### 1. WHAT IS GEM-Cutter?

Gem-cutter is a high-performance bioinformatic library highly optimized for GPUs. The techniques included on the library allow scaling for large problems, using a fine-grain parallelism and advanced programming techniques for GPUs. The library comprises the most used building blocks used for DNA sequence mapping and alignment software. All the algorithms are highly optimized using custom structures and low-level optimizations specifically for each Nvidia GPU architecture (from Fermi to Volta). The library is hiding all the GPU specific programming details using a transparent message passing programming model. The core library is fully integrated and tested with GEM3-mapper, although all its interfaces are well-thought and generalized to be used in other bioinformatics applications. It is self-contained including all the necessary mechanisms for this correct integration on other tools; as dynamic load balance schedulers, I/O modules, memory and compute resource allocators, batch and buffering interfaces, multi-GPU support and wrappers for automatic asynchronous transfers. Gem-cutter is distributed under GPLv3, all the improvements and contributions to this repository are highly appreciated and welcome.

## 2. GETTING STARTED

The Gem-cutter library is fully integrated with GEM3-mapper providing support for multi-GPUs systems, highly optimized specifically for each Nvidia architecture from Fermi to Volta and tested on x86, PPC & aarch64.

A full-fledged last version of the Gem3-GPU mapper can be downloaded from the next repository:
https://github.com/achacond/gem3-mapper

The library can be downloaded as a single module to be integrated into other applications with the next steps:
```
git clone --recursive https://github.com/achacond/gem-cutter.git gem-cutter
cd gem-cutter
make clean all
```

## 3. GPU ACCELERATED MODULES

The core of the Gem-cutter library is a collection of different bioinformatic GPU primitives, CPU-GPU schedulers, resource managers, index and text builders and I/O modules with the aim to make fully transparent to the user the usage of the GPUs, hiding all the system complexity. 
The basic building blocks currently included are described in the next list:

```
* Static Seed FM-index Search
* Adaptative Seed FM-index Search
* Suffix position FM-index decodification
* Candidates position Suffix-Array decodification
* K-mer counting filtering
* Bitparallel Myers edit distance filtering
* Bitparallel Myers edit distance alignment 
* Smith & Waterman gotoh alignment
```

## 4. CONTRIBUTORS

Developers: 
* Alejandro Chacon, <alejandro.chacond@gmail.com> (Main Developer)

Special Contributors:
* Santiago Marco-Sola, <santiagomsola@gmail.com>
* Oriol Mula Valls, <oriol.mula@gmail.com>

Special Thanks:
* Nvidia: Nuno Subtil, Jacopo Pantaleoni, Mark Berguer, Jonathan Cohen
* Thanks to all the friends supporting this project.

## 5. REPORTING BUGS

Feedback and bug reporting it is highly appreciated. 
Please report any issue or suggestion on github or by email (alejandro.chacond@gmail.com) 

## 6. LICENSE AND CITATION

GEM-Cutter is distributed under GPLv3 license and fully free, technical support and bug fixing is freely provided by the authors and the community.

```
Boosting the FM-index on the GPU: effective techniques to mitigate random memory accesses
Alejandro Chacón, Santiago Marco-Sola, Antonio Espinosa, Paolo Ribeca, Juan Carlos Moure

Thread-cooperative, bit-parallel computation of Levenshtein distance on GPU 
Alejandro Chacón, Santiago Marco-Sola, Antonio Espinosa, Paolo Ribeca, Juan Carlos Moure
```


