# reed_solomon

**Collection of FEC code for testing**

1. Reed_solomon code using c implementation

		build : go build -o fec_go fec.go 

  	execution : ./fec_go <options> <inputFile>

  	example : ./fec_go -d 10 -p 3 -b 1024 100m.txt
   		-d : Number of shards to split the data into (default: 10)
   		-p : Number of parity shards (default: 3)
   		-b : Size of one block (default: 1024)
  
2. Simple xor code using c implementation 
  
		Set SIMD value to value other than 1
		--> #define SIMD 0 
  
		build :
  
		execution :
  
  
3. SIMD xor code using c implementation

		Set SIMD value to 1
  	--> #define SIMD 1
  
  	build :
  
  	execution :
