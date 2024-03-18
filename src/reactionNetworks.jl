using Catalyst
# Reactions to consider

closed_competition = @reaction_network begin
    (k1, k1r), A1 <--> A2
    (k2, k2r), A2 <--> A3 
    (k3, k3r), A3 <--> A4
    (k4, k4r), A4 <--> A1 + A1
    (kb1, kb1r), B1 <--> B2 
    (kb2, kb2r), B2 <--> B3 
    (kb3, kb3r), B3 <--> B4 
    (kb4, kb4r), B4 <--> B1 + B1
end

oneStep_open = @reaction_network begin
    (k1, k1r), A1 + 4F <--> A2
    (k2, k2r), A2 <--> A3 
    (k3, k3r), A3 <--> A4
    (k4, k4r), A4 <--> A1 + A1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
end

oneStep_closed = @reaction_network begin
    (k1, k1r), A1 + 4F <--> A2
    (k2, k2r), A2 <--> A3 
    (k3, k3r), A3 <--> A4
    (k4, k4r), A4 <--> A1 + A1
end

twoStep_open = @reaction_network begin
    (k1, k1r), A1 + 2F <--> A2 
    (k2, k2r), A2 + 2F <--> A3 
    (k3, k3r), A3 <--> A4 
    (k4, k4r), A4 <--> A1 + A1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
end

twoStep_closed = @reaction_network begin
    (k1, k1r), A1 + 2F <--> A2 
    (k2, k2r), A2 + 2F <--> A3 
    (k3, k3r), A3 <--> A4 
    (k4, k4r), A4 <--> A1 + A1
end

threeStep_open = @reaction_network begin
    (k1, k1r), A1 + 2F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 <--> A1 + A1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
end

threeStep_closed = @reaction_network begin
    (k1, k1r), A1 + 2F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 <--> A1 + A1
end

autocatalysis_closed = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
end

autocatalysis_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
end

competitive_A4B4_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + F <--> B2
    (b2, b2r), B2 + F <--> B3
    (b3, b3r), B3 + F <--> B4
    (b4, b4r), B4 + F <--> B1 + B1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
    (f3, f3r), B1 <--> ∅
end

competitive_A4B4_closed = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + F <--> B2
    (b2, b2r), B2 + F <--> B3
    (b3, b3r), B3 + F <--> B4
    (b4, b4r), B4 + F <--> B1 + B1
end

competitive_A4B3_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 2F <--> B2
    (b2, b2r), B2 + F <--> B3
    (b3, b3r), B3 + F <--> B4
    (b4, b4r), B4 <--> B1 + B1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
    (f3, f3r), B1 <--> ∅
end

competitive_A4B3_closed = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 2F <--> B2
    (b2, b2r), B2 + F <--> B3
    (b3, b3r), B3 + F <--> B4
    (b4, b4r), B4 <--> B1 + B1
end

competitive_A4B2_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 2F <--> B2
    (b2, b2r), B2 + F <--> B3
    (b3, b3r), B3 + F <--> B4
    (b4, b4r), B4 <--> B1 + B1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
    (f3, f3r), B1 <--> ∅
end

competitive_A4B2_closed = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 2F <--> B2
    (b2, b2r), B2 + 2F <--> B3
    (b3, b3r), B3 <--> B4
    (b4, b4r), B4 <--> B1 + B1
end

competitive_A4B2_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 2F <--> B2
    (b2, b2r), B2 + 2F <--> B3
    (b3, b3r), B3 <--> B4
    (b4, b4r), B4 <--> B1 + B1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
    (f3, f3r), B1 <--> ∅
end

competitive_A4B1_closed = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 4F <--> B2
    (b2, b2r), B2 <--> B3
    (b3, b3r), B3 <--> B4
    (b4, b4r), B4 <--> B1 + B1
end

competitive_A4B1_open = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    (b1, b1r), B1 + 4F <--> B2
    (b2, b2r), B2 <--> B3
    (b3, b3r), B3 <--> B4
    (b4, b4r), B4 <--> B1 + B1
    (f1, f1r), ∅ <--> F
    (f2, f2r), A1 <--> ∅
    (f3, f3r), B1 <--> ∅
end

ac_competition_d = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    k4, A4 + F --> A1 + A1 
    k5, A1 --> 4F
    (kb1, kb1r), B1 + F <--> B2 
    (kb2, kb2r), B2 + F <--> B3 
    (kb3, kb3r), B3 + F <--> B4 
    (kb4, kb4r), B4 + F <--> B5 
    kb5, B5 + F --> B1 + B1
    kb6, B1 --> 5F
end

ac_competition_d2 = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4, k4r), A4 + F <--> A1 + A1 
    k5, A1 --> 4F
    (kb1, kb1r), B1 + F <--> B2 
    (kb2, kb2r), B2 + F <--> B3 
    (kb3, kb3r), B3 + F <--> B4 
    (kb4, kb4r), B4 + F <--> B5 
    (kb5, kb5r), B5 + F <--> B1 + B1
    kb6, B1 --> 5F
end

ac_competition_d3 = @reaction_network begin
    (k1, k1r), A1 + F <--> A2 
    (k2, k2r), A2 + F <--> A3 
    (k3, k3r), A3 + F <--> A4 
    (k4), A4 + F --> A1 + A1 
    (kb1, kb1r), B1 + F <--> B2 
    (kb2, kb2r), B2 + F <--> B3 
    (kb3, kb3r), B3 + F <--> B4 
    (kb4, kb4r), B4 + F <--> B5 
    (kb5), B5 + F --> B1 + B1
end
