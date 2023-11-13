using Catalyst

oneStep_open = @reaction_network begin
    (k1, k1r), 4F + A <--> B
    (k2, k2r), B <--> C 
    (k3, k3r), C <--> D
    (k4, k4r), D <--> A + A
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
end

oneStep_closed = @reaction_network begin
    (k1, k1r), 4F + A <--> B
    (k2, k2r), B <--> C 
    (k3, k3r), C <--> D
    (k4, k4r), D <--> A + A
end

twoStep_open = @reaction_network begin
    (k1, k1r), 2F + A <--> AF2 
    (k2, k2r), 2F + AF2 <--> AA 
    (k3, k3r), AA <--> AA2 
    (k4, k4r), AA2 <--> A + A
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
end

twoStep_closed = @reaction_network begin
    (k1, k1r), 2F + A <--> AF2 
    (k2, k2r), 2F + AF2 <--> AA 
    (k3, k3r), AA <--> AA2 
    (k4, k4r), AA2 <--> A + A
end

threeStep_open = @reaction_network begin
    (k1, k1r), 2F + A <--> AF2 
    (k2, k2r), F + AF2 <--> AF3 
    (k3, k3r), F + AF3 <--> AF4 
    (k4, k4r), AF4 <--> A + A
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
end

threeStep_closed = @reaction_network begin
    (k1, k1r), 2F + A <--> AF2 
    (k2, k2r), F + AF2 <--> AF3 
    (k3, k3r), F + AF3 <--> AF4 
    (k4, k4r), AF4 <--> A + A
end

autocatalysis_closed = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
end

autocatalysis_open = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
end

competitive_A4B3_closed = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 2F <--> BF2
    (b2, b2r), BF2 + F <--> BF3
    (b3, b3r), BF3 + F <--> BF4
    (b4, b4r), BF4 <--> B + B
end

competitive_A4B2_open = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 2F <--> BF2
    (b2, b2r), BF2 + F <--> BF3
    (b3, b3r), BF3 + F <--> BF4
    (b4, b4r), BF4 <--> B + B
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
    (f3, f3r), B <--> ∅
end

competitive_A4B2_closed = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 2F <--> BF2
    (b2, b2r), BF2 + 2F <--> BB
    (b3, b3r), BB <--> BB2
    (b4, b4r), BB2 <--> B + B
end

competitive_A4B2_open = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 2F <--> BF2
    (b2, b2r), BF2 + 2F <--> BB
    (b3, b3r), BB <--> BB2
    (b4, b4r), BB2 <--> B + B
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
    (f3, f3r), B <--> ∅
end

competitive_A4B1_closed = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 4F <--> BB
    (b2, b2r), BB <--> BBa
    (b3, b3r), BBa <--> BBb
    (b4, b4r), BBb <--> B + B
end

competitive_A4B1_open = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    (b1, b1r), B + 4F <--> BBa
    (b2, b2r), BBa <--> BBb
    (b3, b3r), BBb <--> BBc
    (b4, b4r), BBc <--> B + B
    (f1, f1r), ∅ <--> F
    (f2, f2r), A <--> ∅
    (f3, f3r), B <--> ∅
end

ac_competition_d = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    k4, AF3 + F --> A + A 
    k5, A --> 4F
    (kb1, kb1r), B + F <--> BF 
    (kb2, kb2r), BF + F <--> BF2 
    (kb3, kb3r), BF2 + F <--> BF3 
    (kb4, kb4r), BF3 + F <--> BF4 
    kb5, BF4 + F --> B + B
    kb6, B --> 5F
end

ac_competition_d2 = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4, k4r), AF3 + F <--> A + A 
    k5, A --> 4F
    (kb1, kb1r), B + F <--> BF 
    (kb2, kb2r), BF + F <--> BF2 
    (kb3, kb3r), BF2 + F <--> BF3 
    (kb4, kb4r), BF3 + F <--> BF4 
    (kb5, kb5r), BF4 + F <--> B + B
    kb6, B --> 5F
end

ac_competition_d3 = @reaction_network begin
    (k1, k1r), A + F <--> AF 
    (k2, k2r), AF + F <--> AF2 
    (k3, k3r), AF2 + F <--> AF3 
    (k4), AF3 + F --> A + A 
    (kb1, kb1r), B + F <--> BF 
    (kb2, kb2r), BF + F <--> BF2 
    (kb3, kb3r), BF2 + F <--> BF3 
    (kb4, kb4r), BF3 + F <--> BF4 
    (kb5), BF4 + F --> B + B
end

