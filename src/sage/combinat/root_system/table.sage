# A2
branching_rule("A2","A1","levi")

# A3
branching_rule("A3","A2","levi")
branching_rule("A3","A1xA1","tensor")
branching_rule("A3","C2","symmetric")
branching_rule("A3","A1xA1","levi")

# A4
branching_rule("A4","A3","levi")
branching_rule("A4","B2","symmetric")
branching_rule("A4","A1xA2","levi")

# A5
branching_rule("A5","A4","levi")
branching_rule("A5","D3","symmetric")*branching_rule("D3","A3","isomorphic")
branching_rule("A5","A3(0,1,0)","plethysm") # alternative
branching_rule("A5","C3","symmetric")
branching_rule("A5","A2(2,0)","plethysm")
branching_rule("A5","A1xA2","tensor")
branching_rule("A5","A1xA3","levi")
branching_rule("A5","A2xA2","levi")

# A6
branching_rule("A6","A5","levi")
branching_rule("A6","B3","symmetric")
branching_rule("A6","A1xA4","levi")
branching_rule("A6","A2xA3","levi")

# A7
branching_rule("A7","A6","levi")
branching_rule("A7","C4","symmetric")
branching_rule("A7","D4","symmetric")
branching_rule("A7","A1xA3","tensor")
branching_rule("A7","A1xA5","levi")
branching_rule("A7","A2xA4","levi")
branching_rule("A7","A3xA3","levi")

# A8
branching_rule("A8","A7","levi")
branching_rule("A8","B4","symmetric")
branching_rule("A8","A2xA2","tensor")
branching_rule("A8","A1xA6","levi")
branching_rule("A8","A2xA5","levi")
branching_rule("A8","A3xA4","levi")

# B3
branching_rule("B3","G2","miscellaneous")
branching_rule("B3","D3","extended")*branching_rule("D3","A3","isomorphic")
branching_rule("B3","D2xB1","orthogonal_sum")*branching_rule("D2xB1","A1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("B1","A1","isomorphic")])

# B4
branching_rule("B4","D4","extended")
branching_rule("B4","A1","symmetric_power")
branching_rule("B4","B1xB1","tensor")*branching_rule("B1xB1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("B1","A1","isomorphic")])
branching_rule("B4","D2xB2","extended")*branching_rule("D2xB2","A1xA1xB2",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("B4","B1xD3","extended")*branching_rule("B1xD3","A1xA3",[branching_rule("B1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])

# B5
branching_rule("B5","D5","extended")
branching_rule("B5","A1","symmetric_power")
branching_rule("B5","D2xB3","extended")*branching_rule("D2xB3","A1xA2xB3",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("B5","B1xD4","orthogonal_sum")*branching_rule("B1xD4","A1xD4",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("B5","D3xB2","orthogonal_sum")*branching_rule("D3xB2","A3xB2",[branching_rule("D3","A3","isomorphic"),"identity"])

# B6
branching_rule("B6","D6","extended")
branching_rule("B6","A1","symmetric_power")
branching_rule("B6","D2xB4","orthogonal_sum")*branching_rule("D2xB4","A1xA1xB4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("B6","B1xD5","orthogonal_sum")*branching_rule("B1xD5","A1xD5",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("B6","D3xB3","orthogonal_sum")*branching_rule("D3xB3","A3xB3",[branching_rule("D3","A3","isomorphic"),"identity"])
branching_rule("B6","B2xD4","orthogonal_sum")

# B7
branching_rule("B7","D7","extended")
branching_rule("B7","A3(1,0,1)","plethysm")
branching_rule("B7","A1","symmetric_power")
branching_rule("B7","B1xB2","tensor")*branching_rule("B1xB2","A1xB2",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("B7","B1xD6","extended")*branching_rule("B1xD6","A1xD6",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("B7","D2xB5","extended")*branching_rule("D2xB5","A1xA1xB5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("B7","B2xD5","orthogonal_sum")
branching_rule("B7","D3xB4","orthogonal_sum")*branching_rule("D3xB4","A3xB4",[branching_rule("D3","A3","isomorphic"),"identity"])
branching_rule("B7","B3xD4","orthogonal_sum")

# B8
branching_rule("B8","D8","extended")
branching_rule("B8","A1","symmetric_power")
branching_rule("B8","B1xD7","orthogonal_sum")*branching_rule("B1xD7","A1xD7",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("B8","D2xB6","orthogonal_sum")*branching_rule("D2xB6","A1xA1xB6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("B8","B2xD6","orthogonal_sum")
branching_rule("B8","D3xB5","orthogonal_sum")*branching_rule("D3xB5","A3xB5",[branching_rule("D3","A3","isomorphic"),"identity"])
branching_rule("B8","B3xD5","orthogonal_sum")
branching_rule("B8","B4xD4","orthogonal_sum")

# C2
branching_rule("C2","A1","symmetric_power")
branching_rule("C2","C1xC1","orthogonal_sum")*branching_rule("C1xC1","A1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])

# C3
branching_rule("C3","A2","levi")
branching_rule("C3","A1","symmetric_power")
branching_rule("C3","B1xC1","tensor")*branching_rule("B1xC1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])
branching_rule("C3","C1xC2","orthogonal_sum")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])

# C4
branching_rule("C4","A3","levi")
branching_rule("C4","A1","symmetric_power")
branching_rule("C4","C1xC3","orthogonal_sum")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C4","C2xC2","orthogonal_sum")
branching_rule("C4","C1xD2","tensor")*branching_rule("C1xD2","A1xA1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])

# C5
branching_rule("C5","A4","levi")
branching_rule("C5","A1","symmetric_power")
branching_rule("C5","C1xC4","orthogonal_sum")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C5","C2xC3","orthogonal_sum")
branching_rule("C5","C1xB2","tensor")*branching_rule("C1xB2","A1xB2",[branching_rule("C1","A1","isomorphic"),"identity"])

# C6
branching_rule("C6","A5","levi")
branching_rule("C6","A1","symmetric_power")
branching_rule("C6","C1xD3","tensor")*branching_rule("C1xD3","A1xA3",[branching_rule("C1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])
branching_rule("C6","B1xC2","tensor")*branching_rule("B1xC2","A1xC2",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("C6","C1xC5","orthogonal_sum")*branching_rule("C1xC5","A1xC5",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C6","C2xC4","orthogonal_sum")
branching_rule("C6","C3xC3","orthogonal_sum")

# C7
branching_rule("C7","A6","levi")
branching_rule("C7","A1","symmetric_power")
branching_rule("C7","C1xB3","tensor")*branching_rule("C1xB3","A1xB3",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C7","C1xC6","orthogonal_sum")*branching_rule("C1xC6","A1xC6",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C7","C2xC5","orthogonal_sum")
branching_rule("C7","C3xC4","orthogonal_sum")
branching_rule("C7","C3(0,0,1)","plethysm") # overlooked by Patera and McKay

# C8
branching_rule("C8","A7","levi")
branching_rule("C8","A1","symmetric_power")
branching_rule("C8","C2(1,1)","plethysm")
branching_rule("C8","C1xD4","tensor")*branching_rule("C1xD4","A1xD4",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C8","C1xC7","orthogonal_sum")*branching_rule("C1xC7","A1xC7",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("C8","C2xC6","orthogonal_sum")
branching_rule("C8","C3xC5","orthogonal_sum")
branching_rule("C8","C4xC4","orthogonal_sum")

# D4
branching_rule("D4","B3","symmetric")
branching_rule("D4","A2(1,1)","plethysm")
branching_rule("D4","C1xC2","tensor")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("D4","D2xD2","orthogonal_sum")*branching_rule("D2xD2","A1xA1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])

# D5
branching_rule("D5","A4","levi")
branching_rule("D5","B4","symmetric")
branching_rule("D5","C2(2,0)","plethysm")
branching_rule("D5","D2xD3","orthogonal_sum")*branching_rule("D2xD3","A1xA1xA3",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D3","A3","isomorphic")])
branching_rule("D5","B1xB3","orthogonal_sum")*branching_rule("B1xB3","A1xA3",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("D5","B2xB2","orthogonal_sum")

# D6
branching_rule("D6","A5","levi")
branching_rule("D6","B5","symmetric")
branching_rule("D6","C1xC3","tensor")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("D6","D2xD4","orthogonal_sum")*branching_rule("D2xD4","A1xA1xD4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("D6","D3xD3","orthogonal_sum")*branching_rule("D3xD3","A3xA3",[branching_rule("D3","A3","isomorphic"),branching_rule("D3","A3","isomorphic")])
branching_rule("D6","B1xB4","orthogonal_sum")*branching_rule("B1xB4","A1xB4",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("D6","B2xB3","orthogonal_sum")
branching_rule("D6","B1xD2","tensor")*branching_rule("B1xD2","A1xA1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])

# D7
branching_rule("D7","A6","levi")
branching_rule("D7","B6","symmetric")
branching_rule("D7","C3(0,1,0)","plethysm")
branching_rule("D7","C2(0,2)","plethysm")
branching_rule("D7","G2(0,1)","plethysm")
branching_rule("D7","D2xD5","orthogonal_sum")*branching_rule("D2xD5","A1xA1xD5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("D7","D3xD4","orthogonal_sum")*branching_rule("D3xD4","A3xD4",[branching_rule("D3","A3","isomorphic"),"identity"])
branching_rule("D7","B1xB5","orthogonal_sum")*branching_rule("B1xB5","A1xB5",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("D7","B2xB4","orthogonal_sum")
branching_rule("D7","B3xB3","orthogonal_sum")

# D8
branching_rule("D8","A7","levi")
branching_rule("D8","B7","symmetric")
branching_rule("D8","B4(0,0,0,1)","plethysm")
branching_rule("D8","C1xC4","tensor")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])
branching_rule("D8","D2xD6","orthogonal_sum")*branching_rule("D2xD6","A1xA1xD6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])
branching_rule("D8","D3xD5","orthogonal_sum")*branching_rule("D3xD5","A3xD5",[branching_rule("D3","A3","isomorphic"),"identity"])
branching_rule("D8","D4xD4","orthogonal_sum")
branching_rule("D8","B1xB6","orthogonal_sum")*branching_rule("B1xB6","A1xB6",[branching_rule("B1","A1","isomorphic"),"identity"])
branching_rule("D8","B2xB5","orthogonal_sum")
branching_rule("D8","B3xB4","orthogonal_sum")
branching_rule("D8","C2xC2","tensor")

# G2
branching_rule("G2","A2","extended")
branching_rule("G2","A1","i")
branching_rule("G2","A1xA1","extended")

# F4
branching_rule("F4","B4","extended")
branching_rule("F4","A1","ii")
branching_rule("F4","A1xG2","miscellaneous")
branching_rule("F4","A1xC3","extended")
branching_rule("F4","A2xA2","extended")

# E6
branching_rule("E6","D5","levi")
branching_rule("E6","C4","symmetric")
branching_rule("E6","F4","symmetric")
branching_rule("E6","A2","miscellaneous")
branching_rule("E6","G2","miscellaneous")
branching_rule("E6","A2xG2","miscellaneous")
branching_rule("E6","A1xA5","extended")
branching_rule("E6","A2xA2xA2","extended")

# E7
branching_rule("E7","A7","extended")
branching_rule("E7","E6","levi")
branching_rule("E7","A2","miscellaneous")
branching_rule("E7","A1","iii")
branching_rule("E7","A1","iv")
branching_rule("E7","A1xF4","miscellaneous")
branching_rule("E7","G2xC3","miscellaneous")
branching_rule("E7","A1xG2","miscellaneous")
branching_rule("E7","A1xA1","miscellaneous")
branching_rule("E7","A1xD6","extended")
branching_rule("E7","A5xA2","extended")

# E8
branching_rule("E8","A4xA4","extended")
branching_rule("E8","G2xF4","miscellaneous")
branching_rule("E8","E6xA2","extended")
branching_rule("E8","E7xA1","extended")
branching_rule("E8","D8","extended")
branching_rule("E8","A8","extended")
branching_rule("E8","B2","miscellaneous")
branching_rule("E8","A1xA2","miscellaneous")
branching_rule("E8","A1","v")
branching_rule("E8","A1","vi")
branching_rule("E8","A1","vii")

