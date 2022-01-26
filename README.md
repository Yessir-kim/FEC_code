# FEC

**Collection of FEC code for testing**

1. Reed_solomon code using c implementation

		build : go build -o fec_go fec_rs.go 

  		execution : ./fec_go <options> <inputFile>

  		example : ./fec_go -d 10 -p 3 -b 1024 100m.txt
   			-d : Number of shards to split the data into (default: 10)
   			-p : Number of parity shards (default: 3)
   			-b : Size of one block (default: 1024)
  
2. Simple xor code using c implementation 
  
		Set SIMD value to value other than 1
		--> #define SIMD 0 in fec_xor.go 
  
		build : go build -o fec_go fec_xor.go
  
		execution : ./fec_go <options> <inputFile>

		example : ./fec_go -r 32 -c 1024 100m.txt
			-r : Size of row
			-c : Size of column
			--> block size = row * col
  
  
3. SIMD xor code using c implementation

		Set SIMD value to 1
  		--> #define SIMD 1 in fec_xor.go
  
  		build : go build -o fec_go fec_xor.go
  
  		execution : ./fec_go <options> <inputFile>

		example : ./fec_go -r 32 -c 1024 100m.txt
			-r : Size of row
			-c : Size of column
			--> block size = row * col
