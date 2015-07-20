import all
import grs

linear_code.AbstractLinearCode._registered_encoders["GeneratorMatrix"] = linear_code.LinearCodeGeneratorMatrixEncoder
linear_code.AbstractLinearCode._registered_decoders["Syndrome"] = linear_code.LinearCodeSyndromeDecoder
linear_code.AbstractLinearCode._registered_decoders["NearestNeighbor"] = linear_code.LinearCodeNearestNeighborDecoder

grs.GeneralizedReedSolomonCode._registered_encoders["EvaluationVector"] = grs.GRSEvaluationVectorEncoder
grs.GeneralizedReedSolomonCode._registered_encoders["EvaluationPolynomial"] = grs.GRSEvaluationPolynomialEncoder
grs.GeneralizedReedSolomonCode._registered_decoders["Syndrome"] = linear_code.LinearCodeSyndromeDecoder
grs.GeneralizedReedSolomonCode._registered_decoders["NearestNeighbor"] = linear_code.LinearCodeNearestNeighborDecoder
