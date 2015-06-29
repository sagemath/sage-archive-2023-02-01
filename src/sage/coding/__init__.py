import all
linear_code.AbstractLinearCode._registered_decoders["Syndrome"] = linear_code.LinearCodeSyndromeDecoder
linear_code.AbstractLinearCode._registered_decoders["NearestNeighbor"] = linear_code.LinearCodeNearestNeighborDecoder
