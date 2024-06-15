# ------------------------------------------------------------------------------
# ---- 1. INTRODUCTION AND CONTENTS --------------------------------------------

# Packages used in this code.
using Dates
using IntervalArithmetic    # 0.22.14 
using BenchmarkTools        # 1.5.0 

function dateAndTime(string) # The string is the message shown while printing the human-readable date and time

    # time() outputs the number of seconds since the Unix epoch (January 1, 1970) so can easily be used for calculations.
    currentTime = time()

    # now() outputs a human-readable date and time (e.g. "2022-04-16T03:38:17.847") so we print it.
    currentDate = now()
    println("$string: $currentDate")

    return [currentTime currentDate]

end

# Converts a floating-point number to an IntervalArithmetic "number." 
I(x) = @interval(x) 

LO(I) = IntervalArithmetic.inf(I)
HI(I) = IntervalArithmetic.sup(I)

# Integer powers in IntervalArithmetic are much faster with repeated multiplication than using "^"; the latter converts the numbers into BigFloat's. 
function P(X, power) # power must be a POSITIVE INTEGER 
    result = I(1) 
    for i = 1 : power 
        result = result * I(X) 
    end 
    return result 
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 2. FORMULAS FOR THE 24-CELL ---------------------------------------------

# BASIC VOLUME AND NUMBER OF POINTS FORMULAS 

vol_sausage(n) = 2 * (n - 1) * (4 / 3 * π) + 1 / 2 * π^2 
        g_Y(m) = 4 * m^4 + 8 * m^3 + 8 * m^2 + 4 * m + 1 
      vol_Y(m) = 32 * m^4 

Vol_sausage(N) = I(2) * (N - I(1)) * (I(4) / I(3) * I(π)) + I(1) / I(2) * I(π)*I(π) 
        G_Y(M) = I(4) * P(M, 4) + I(8) * P(M, 3) + I(8) * P(M, 2) + I(4 * M + 1) 
      Vol_Y(M) = I(32) * P(M, 4) 



# SUMS OF POWERS OF TRUNCATIONS 

h_(h) = [h[1] + h[2] + h[3]
         h[1]^2 + h[2]^2 + h[3]^2
         h[1]^3 + h[2]^3 + h[3]^3
         h[1]^4 + h[2]^4 + h[3]^4] 

H_(H) = [H[1] + H[2] + H[3] 
         P(H[1], 2) + P(H[2], 2) + P(H[3], 2) 
         P(H[1], 3) + P(H[2], 3) + P(H[3], 3) 
         P(H[1], 4) + P(H[2], 4) + P(H[3], 4) ]



# THE NUMBER OF POINTS OF A FACET-TRUNCATED 24-CELL 

function g_t3h_Y(m, h)

    # Components of g_t3h_Y. 
    c = 8 - 2 * h_(h)[1] / 3
    f = 8 - (h_(h)[2] + h_(h)[1])
    e = 4 - (2 * h_(h)[3] + 3 * h_(h)[2] + 2 * h_(h)[1]) / 3
    v = 1 - (-2 * h_(h)[4] + 2 * h_(h)[3] + 5 * h_(h)[2] + h_(h)[1]) / 6 

    # Number of points of t3h_Y. 
    return 4 * m^4 + c * m^3 + f * m^2 + e * m + v 

end 

function G_t3h_Y(M, H)

    # Components of G_t3h_Y. 
    C = I(8) - I(2) * H_(H)[1] / I(3) 
    F = I(8) - (H_(H)[2] + H_(H)[1]) 
    E = I(4) - (I(2) * H_(H)[3] + I(3) * H_(H)[2] + I(2) * H_(H)[1]) / I(3) 
    V = I(1) - (-I(2) * H_(H)[4] + I(2) * H_(H)[3] + I(5) * H_(H)[2] + H_(H)[1]) / I(6) 

    # Number of points of t3h_Y. 
    return I(4) * P(M, 4) + C * P(M, 3) + F * P(M, 2) + E * M + V 

end 



# VOLUME OF THE FACET-TRUNCATED 24-CELL 

function vol_t3h_Y(m, h) 

    # Components of vol(t_h_3(Y_m) + B^4). 
    c = 64 * √2 - 16 / 3 * h_(h)[1] 
    f = 16 * √3 * π - 8 * √2 * h_(h)[1] - 8 * h_(h)[2] 
    e = 64 * (3 * acos(1 / 3) - π) - 4 * √3 * π / 3 * h_(h)[1] - 8 * √2 * h_(h)[2] - 16 / 3 * h_(h)[3] 
    v = 1 / 2 * π^2 + (64 * atan(3 - 2 * √2) - 24 * (3 * acos(1 / 3) - π)) * h_(h)[1] + (18 - 14 * √3) * π / 3 * h_(h)[2] - 8 * √2 / 3 * h_(h)[3] + 8 / 3 * h_(h)[4] 

    # Volume of t_h_3(Y_m) + B^4. 
    return vol_Y(m) + c * m^3 + f * m^2 + e * m + v 

end 

function Vol_t3h_Y(M, H) 

    # Components of vol(t_h_3(Y_m) + B^4). 
    C = I(64) * √I(2) - I(16) / I(3) * H_(H)[1] 
    F = I(16) * √I(3) * I(π) - I(8) * √I(2) * H_(H)[1] - I(8) * H_(H)[2] 
    E = I(64) * (I(3) * acos(I(1) / I(3)) - I(π)) - I(4) * √I(3) * I(π) / I(3) * H_(H)[1] - I(8) * √I(2) * H_(H)[2] - I(16) / I(3) * H_(H)[3] 
    V = I(1) / I(2) * P(π, 2) + (I(64) * atan(I(3) - I(2) * √I(2)) - I(24) * (I(3) * acos(I(1) / I(3)) - I(π))) * H_(H)[1] + (I(18) - I(14) * √I(3)) * I(π) / I(3) * H_(H)[2] - I(8) * √I(2) / I(3) * H_(H)[3] + I(8) / I(3) * H_(H)[4] 

    # Volume of t_h_3(Y_m) + B^4. 
    return Vol_Y(M) + C * P(M, 3) + F * P(M, 2) + E * M + V 

end 



# The following functions are "for show" only, so do not need to use interval arithmetic. 

# The density of the four-dimensional sausage with n points. 
δ4_sausage(n) = (n * (1 / 2 * π^2)) / vol_sausage(n) 

# The density of the truncated 24-cell. 
δ4_t3h_Y(m, h) = (g_t3h_Y(m, h) * (1 / 2 * π^2)) / vol_t3h_Y(m, h) 

# The "approximate" density of the truncated 24-cell with k points removed. 
δ̃4_t3h_Y(m, h, k) = ((g_t3h_Y(m, h) - k) * (1 / 2 * π^2)) / vol_t3h_Y(m, h) 



# The maximum number of points that can be removed from a packing to keep it denser than the sausage. 
ṽ(m, h) = vol_sausage(g_t3h_Y(m, h)) - vol_t3h_Y(m, h) 
r̃(m, h) = ceil(ṽ(m, h) / (2 * (4 / 3 * π)) - 1) 
Ṽ(M, H) = Vol_sausage(G_t3h_Y(M, H)) - Vol_t3h_Y(M, H) 
R̃(M, H) = ceil(Ṽ(M, H) / (I(2) * (I(4) / I(3) * I(π))) - I(1)) 

# Let n = p_k(t_h_3(Z_m)). This function outputs a range of n for which p_k(t_h_3(Z_m)) is denser than S_n. 
# If p_k(t_h_3(Z_m)) is LESS dense than S_n then it will, of course, be incorrect. 
L̃(m, h) = [g_t3h_Y(m, h) - r̃(m, h)   g_t3h_Y(m, h)]
L̃_IA(M, H) = [G_t3h_Y(M, H) - R̃(M, H)   G_t3h_Y(M, H)]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ---- 3. DENSITY COMPARISONS --------------------------------------------------

# We are not going to list all dimensions for which 
# δ̃(p_k(t_h_3(Z_m))) > δ((S_n)^d). 
# Instead, we display ranges for which that inequality is true and false. 
# That will make the ranges more reasonable. 

# This function will take an input of n and determine if some packing in Z3 with n spheres is δ̃-denser than the sausage. 
function check_n(n, details = false, fast = false) 
    
    if details == true 
        println("")
        startDate = dateAndTime("         Start date")
        startDate
        println("")
    end 
    is_denser_than_sausage = 0 
    h_max(m) = floor((m - 1) / 2)

    min_m = 1 
    max_m = 1 

    if fast == true # Shortcuts for when fast == true 

        while g_Y(min_m) < n 
            min_m += 1 
        end 

        max_m = 9999

    else # if fast == false 

        while max( g_Y(min_m), (g_Y(min_m) - r̃(min_m, [h_max(min_m) h_max(min_m) h_max(min_m)])) ) < n 
            min_m += 1 
        end 

        while g_t3h_Y(max_m, [h_max(max_m) h_max(max_m) h_max(max_m)]) ≤ n 
            max_m += 1 
        end # Sometimes (for small n) max_m will be less than min_m — that's fine since we only use this code for n ≳ 500,000 anyway 
        max_m -= 1 

    end 

    m = 1 
    loopCounter = 0

    for m = min_m : max_m 
        if details == true 
            print("$m ")
        end 
        for h3 = 0 : h_max(m) 
            for h2 = 0 : h3 
                for h1 = 0 : h2 
                    range_of_packing = L̃(m, [h1 h2 h3])
                    loopCounter += 1 
                    if ( (range_of_packing[1] - 0.5) ≤ n) & (n ≤ (range_of_packing[2] + 0.5) ) == true # NOTE: These inequalities cannot both be true if p_k(t_h_3(Z_m)) is LESS dense than S_n 

                        if details == true 
                            println("")
                            print("$m, [$(Int(h1)) $(Int(h2)) $(Int(h3))] ")
                        end 
                        is_denser_than_sausage += 1 

                        if (is_denser_than_sausage > 0) & (fast == true) 
                            if details == true 
                                println("")
                                println("")
                                endDate = dateAndTime("       End date")
                                endDate
                                println("     Time taken: $(endDate[1] - startDate[1]) seconds")
                                println("Number of loops: $loopCounter")
                                println("          Speed: $(loopCounter / (endDate[1] - startDate[1])) calculations of L per second.")
                                println("")
                            end 
                        
                            return is_denser_than_sausage
                        end 

                    end 
                    
                end 
            end 
        end 
        m += 1 
    end 

    if details == true 
        println("")
        println("")
        endDate = dateAndTime("           End date")
        endDate
        println("     Time taken: $(endDate[1] - startDate[1]) seconds")
        println("Number of loops: $loopCounter")
        println("          Speed: $(loopCounter / (endDate[1] - startDate[1])) calculations of L̃ per second.")
        println("")
    end 

    return is_denser_than_sausage

end 

function Check_n(n, details = false, fast = false) 

    N = I(n) 
    
    if details == true 
        println("")
        startDate = dateAndTime("         Start date")
        startDate
        println("")
    end 
    is_denser_than_sausage = 0 
    h_max(m) = floor((m - 1) / 2)

    min_m = 1 
    max_m = 1 

    if fast == true # Shortcuts for when fast == true 

        while HI(G_Y(min_m)) < LO(N)
            min_m += 1 
        end 

        max_m = 9999

    else # if fast == false 

        while max( HI(G_Y(I(min_m))), HI(G_Y(I(min_m)) - R̃(I(min_m), [h_max(I(min_m)) h_max(I(min_m)) h_max(I(min_m))])) ) < LO(N)
            min_m += 1 
        end 

        while HI(G_t3h_Y(I(max_m), [h_max(I(max_m)) h_max(I(max_m)) h_max(I(max_m))])) ≤ LO(N)
            max_m += 1 
        end # Sometimes (for small n) max_m will be less than min_m — that's fine since we only use this code for n ≳ 500,000 anyway 
        max_m -= 1 

    end 

    m = 1 
    loopCounter = 0

    for m = min_m : max_m 

        M = I(m) 

        if details == true 
            print("$m ")
        end 
        for h3 = 0 : h_max(m) 
            for h2 = 0 : h3 
                for h1 = 0 : h2 
                    range_of_packing = L̃_IA(M, I.([h1 h2 h3]))
                    loopCounter += 1 
                    if ( HI(range_of_packing[1] - I(0.5)) ≤ LO(N)) & (HI(N) ≤ LO(range_of_packing[2] + I(0.5)) ) == true # NOTE: These inequalities cannot both be true if p_k(t_h_3(Y_m)) is LESS dense than S_n 

                        if details == true 
                            println("")
                            print("$m, [$(Int(h1)) $(Int(h2)) $(Int(h3))] ")
                        end 
                        is_denser_than_sausage += 1 

                        if (is_denser_than_sausage > 0) & (fast == true) 
                            if details == true 
                                println("")
                                println("")
                                endDate = dateAndTime("       End date")
                                endDate
                                println("     Time taken: $(endDate[1] - startDate[1]) seconds")
                                println("Number of loops: $loopCounter")
                                println("          Speed: $(loopCounter / (endDate[1] - startDate[1])) calculations of L per second.")
                                println("")
                            end 
                        
                            return is_denser_than_sausage
                        end 

                    end 
                    
                end 
            end 
        end 
        m += 1 
    end 

    if details == true 
        println("")
        println("")
        endDate = dateAndTime("           End date")
        endDate
        println("     Time taken: $(endDate[1] - startDate[1]) seconds")
        println("Number of loops: $loopCounter")
        println("          Speed: $(loopCounter / (endDate[1] - startDate[1])) calculations of L̃ per second.")
        println("")
    end 

    return is_denser_than_sausage

end 

#= 
FFS Julia: 

julia> L̃(19, [0 1 7])[1]
514649.99999999994

julia> L̃(19, [0 1 7])[2]
516835.99999999994
=#

#= Justification for the condition "if ( LO(range_of_packing[1]) ≤ n ) & HI(n ≤ range_of_packing[2]) == true": 
    The calculation of L will output an interval, e.g. 

    L̃(17, [1 3 4])[1] = [338195.99999999994, 338196.00000000006]

    (use L̃(@interval(17), [@interval(1) @interval(3) @interval(4)])[1].lo). 
    The true value is n = 338,196 and we want the claim L̃(17, [1 3 4])[1] ≤ n. 
    Since L̃(17, [1 3 4])[1] outputs, in general, [338,196 - ε, 338,196 + ε], we know that if the upper bound is less than 338,197 then L̃(17, [1 3 4])[1] is definitely ≤ 338,196. 
    So we use the inequality 

    HI(L̃(17, [1 3 4])[1]) < n + 0.5 

    for the interval condition. Similarly, n ≤ L̃(17, [1 3 4])[2] becomes 
    
    n - 0.5 < LO(L̃(17, [1 3 4])[2]). 

=# 



# This function will check all natural numbers from 5 to some number N to see which definitely have packings denser than the sausage and which might not. 
# The output will be the number of spheres and expressed in ranges. 

function check_n_range(n_start, n_end, details = false, fast = false) 
    #= 
    fast MUST be set to FALSE for n_end < 516,946 or the function will NOT terminate! 
    fast should be set to TRUE for n_start ≥ 516,946 to speed up computations by roughly a factor of 8. 
    As a corollary, DO NOT USE any case of n_start < 516,946 ≤ n_end. (n_start = 516,946 is fine.) 
    =#

    println("")
    println("Checking all natural numbers from $n_start to $n_end (inclusive)")
    println("")
    startDate = dateAndTime("         Start date")
    startDate
    currentTime = startDate[1]
    println("")
    if details == true 
        if n_start ≥ 516_946
            if check_n(n_start, false, true) == 0 
                print("        Sausage?: {$n_start, ..., ")
            end 
        else # if n_start < 516_946
            if check_n(n_start, false, false) == 0 
                print("        Sausage?: {$n_start, ..., ")
            end 
        end 
    end 
    last_n_sign = 0 
    this_n_sign = 0 

    for n = n_start : n_end 
        this_n_sign = sign(check_n(n, false, fast)) # If check_n > 0 then set this_n_sign == 1; if check_n_sign = 0 then set this_n_sign == 0 
        if last_n_sign ≠ this_n_sign 
            if details == true 
                if this_n_sign == 1 
                    # If n = n_start then there is no existing "{n, ...," string in the output and you will get a lone "n - 1}" line. Therefore only print a "n - 1}" if there is already a line of the form "{n, ...," in the output. 
                    if n > n_start 
                        println("$(n - 1)}")
                    end 
                    print("Full-dimensional: {$n, ..., ")
                else # if this_n_sign == 0 
                    println("$(n - 1)}")
                    print("        Sausage?: {$n, ..., ")
                end 
            else # if details == false 
                if this_n_sign == 1 
                    print("{$n, ..., ")
                else # if this_n_sign == 0 
                    println("$(n - 1)}")
                end 
            end 
        end 
        last_n_sign = this_n_sign 
        if n % 1_000_000 == 0
            # The speed calculation below assumes that check_n has checked the past 1,000,000 numbers, otherwise it will give incorrect numbers. If that's not the case, then we print only the value of n and not the speed. 
            if currentTime == startDate[1] 
                println("[ n = $n ]") 
            else #currentTime ≠ startDate[1]
                println("[ n = $n    Speed: $(round(1_000_000 / (time() - currentTime), digits = 1)) ]")
            end 
            currentTime = time()
        end 
    end 
    
    if details == true 
        print("$n_end}")
        println("")
    else # details == false 
        if this_n_sign == 1 
            print("$n_end}") # If the sausage "is optimal" then the line doesn't get printed 
            println("")
        end 
    end 

    println("")
    println("The last number checked: n = $n_end")

    println("")
    println("")
    endDate = dateAndTime("           End date")
    endDate
    println("         Time taken: $(round(endDate[1] - startDate[1], digits = 3)) seconds")
    println("Values of n checked: $(n_end - n_start + 1)")
    println("              Speed: $(round((n_end - n_start + 1) / (endDate[1] - startDate[1]), digits = 1)) values of n checked per second.")
    println("")

end 

function Check_n_range(n_start, n_end, details = false, fast = false) 
    #= 
    fast MUST be set to FALSE for n_end < 516,946 or the function will NOT terminate! 
    fast should be set to TRUE for n_start ≥ 516,946 to speed up computations by roughly a factor of 8. 
    As a corollary, DO NOT USE any case of n_start < 516,946 ≤ n_end. (n_start = 516,946 is fine.) 
    =#

    println("")
    println("Checking all natural numbers from $n_start to $n_end (inclusive)")
    println("")
    startDate = dateAndTime("         Start date")
    startDate
    currentTime = startDate[1]
    println("")
    if details == true 
        if n_start ≥ 516_946 
            if Check_n(n_start, false, true) == 0 
                print("        Sausage?: {$n_start, ..., ")
            end 
        else # if n_start < 516_946
            if Check_n(n_start, false, false) == 0 
                print("        Sausage?: {$n_start, ..., ")
            end 
        end 
    end 
    last_n_sign = 0 
    this_n_sign = 0 

    for n = n_start : n_end 
        
        this_n_sign = sign(Check_n(n, false, fast)) 
        # If Check_n > 0 then set this_n_sign == 1; if Check_n = 0 then set this_n_sign == 0. 
        # Check_n outputs either 0 or a positive number, so the above comparison is valid. 
        if last_n_sign ≠ this_n_sign 
            if details == true 
                if this_n_sign == 1 
                    # If n = n_start then there is no existing "{n, ...," string in the output and you will get a lone "n - 1}" line. Therefore only print a "n - 1}" if there is already a line of the form "{n, ...," in the output. 
                    if n > n_start 
                        println("$(n - 1)}")
                    end 
                    print("Full-dimensional: {$n, ..., ")
                else # if this_n_sign == 0 
                    println("$(n - 1)}")
                    print("        Sausage?: {$n, ..., ")
                end 
            else # if details == false 
                if this_n_sign == 1 
                    print("Full-dimensional: {$n, ..., ")
                else # if this_n == 0 
                    println("$(n - 1)}")
                end 
            end 
        end 
        last_n_sign = this_n_sign 
        if n % 1_000_000 == 0
            # The speed calculation below assumes that check_N has checked the past 1,000,000 numbers, otherwise it will give incorrect numbers. If that's not the case, then we print only the value of N and not the speed. 
            if currentTime == startDate[1] 
                println("[ n = $n ]") 
            else #currentTime ≠ startDate[1]
                println("[ n = $n    Speed: $(round(1_000_000 / (time() - currentTime), digits = 1)) ]")
            end 
            currentTime = time()
        end 
    end 
    
    if details == true 
        print("$n_end}")
        println("")
    else # details == false 
        if this_n_sign == 1 
            print("$n_end}") # If the sausage "is optimal" then the line doesn't get printed 
            println("")
        end 
    end 

    println("")
    println("The last number checked: n = $n_end")

    println("")
    println("")
    endDate = dateAndTime("           End date")
    endDate
    println("         Time taken: $(round(endDate[1] - startDate[1], digits = 3)) seconds")
    println("Values of n checked: $(n_end - n_start + 1)")
    println("              Speed: $(round((n_end - n_start + 1) / (endDate[1] - startDate[1]), digits = 1)) values of n checked per second.")
    println("")

end 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------