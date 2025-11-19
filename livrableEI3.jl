using Random

# --------------------------------------------------------------------------- #
#   1. 贪婪修复函数 (Greedy Repair) - [核心改进]
#      逻辑参考 GA_SPP_NEW.jl 
#      
#      不再随机删除，而是先把解中所有的变量按价值(C)从大到小排序。
#      然后像拼图一样，优先把贵的变量放进新解里。
# --------------------------------------------------------------------------- #

"""
    greedy_repair(x_input, A, C)

    接收一个可能不可行(有冲突)的解 x_input。
    返回一个绝对可行且尽可能高价值的解 x_new。
"""
function greedy_repair(x_input, A, C)
    m, n = size(A)
    
    # 1. 找出输入解中所有被选中的变量 (即 x[j] == 1 的 j)
    selected_indices = Int[]
    for j = 1:n
        if x_input[j] == 1
            push!(selected_indices, j)
        end
    end
    
    # 2. 关键步骤：按价值 C[j] 从大到小排序
    #    这样我们优先保留最赚钱的变量
    sort!(selected_indices, by = j -> C[j], rev=true)
    
    # 3. 重建一个全新的、合法的解
    x_new = zeros(Int, n)
    covered_rows = zeros(Bool, m) # 记录哪些行已经被占用了
    
    # 遍历排序后的变量
    for j in selected_indices
        
        # 检查变量 j 是否能放入 (是否与 covered_rows 冲突)
        is_safe = true
        for i = 1:m
            if A[i, j] == 1 && covered_rows[i] == true
                is_safe = false # 冲突了，这一行已经被之前的(更贵的)变量占了
                break
            end
        end
        
        # 如果安全，就放入新解
        if is_safe
            x_new[j] = 1
            # 更新占用的行
            for i = 1:m
                if A[i, j] == 1
                    covered_rows[i] = true
                end
            end
        end
    end
    
    return x_new
end


function crossover_standard(p1, p2)
    n = length(p1)
    point = rand(2:n-1)
    child = zeros(Int, n)
    
    # 前半段用 p1，后半段用 p2
    for j = 1:point
        child[j] = p1[j]
    end
    for j = point+1:n
        child[j] = p2[j]
    end
    return child
end


function crossover_XOR(p1, p2)
    n = length(p1)
    child = zeros(Int, n)
    for j = 1:n
        # Julia 的 xor 函数或 != 都可以
        if p1[j] != p2[j]
            child[j] = 1
        else
            child[j] = 0
        end
    end
    return child
end


function crossover_AND(p1, p2)
    n = length(p1)
    child = zeros(Int, n)
    for j = 1:n
        if p1[j] == 1 && p2[j] == 1
            child[j] = 1
        else
            child[j] = 0
        end
    end
    return child
end


function genetic_algorithm_SPP(C, A, max_time=10.0)
    m, n = size(A)
    start_time = time()
    
    # --- 参数设置 ---
    pop_size = 30      # 种群大小
    p_mutation = 0.1   # 变异概率
    max_generations = 2000 
    
    # --------------------------------
    # 1. 初始化种群
    # --------------------------------
    population = [] 
    
    println("  [GA] 初始化种群 (使用 EI1 贪婪算法)...")
    for k = 1:pop_size
        # 随机 alpha 生成多样化初始解
        rand_alpha = 0.5 + rand() * 0.5 
        x_init, z_init = randomized_greedy_construction(C, A, rand_alpha)
        push!(population, x_init)
    end
    
    # 记录最好的
    best_x = copy(population[1])
    best_z = dot(C, best_x)
    
    # --------------------------------
    # 2. 世代循环
    # --------------------------------
    for gen = 1:max_generations
        
        # 时间检查
        if time() - start_time > max_time
            println("  [GA] 时间耗尽，停止于第 ", gen, " 代。")
            break
        end
        
        new_population = []
        
        # 生成下一代
        while length(new_population) < pop_size
            
            # --- A. 选择 (锦标赛) ---
            # 选爸爸
            i1, i2 = rand(1:pop_size), rand(1:pop_size)
            p1 = dot(C, population[i1]) > dot(C, population[i2]) ? population[i1] : population[i2]
            # 选妈妈
            i3, i4 = rand(1:pop_size), rand(1:pop_size)
            p2 = dot(C, population[i3]) > dot(C, population[i4]) ? population[i3] : population[i4]
            
            # --- B. 交叉 (混合策略) ---
            r_op = rand()
            child_raw = zeros(Int, n)
            
            if r_op < 0.5
                # 50% 概率：普通单点交叉
                child_raw = crossover_standard(p1, p2)
            elseif r_op < 0.75
                # 25% 概率：XOR (异或)
                child_raw = crossover_XOR(p1, p2)
            else
                # 25% 概率：AND (与)
                child_raw = crossover_AND(p1, p2)
            end
            
            # --- C. 变异 (Mutation) ---
            if rand() < p_mutation
                idx = rand(1:n)
                child_raw[idx] = 1 - child_raw[idx] # 翻转 0<->1
            end
            
            # --- D. 贪婪修复 (Greedy Repair) ---
            # 这是关键：修复并优化不可行的子代
            child_feasible = greedy_repair(child_raw, A, C)
            
            # --- E. 局部搜索 (Memetic 强化) ---
            # 调用 EI1 的局部搜索让孩子更强
            child_improved, z_child = local_search(child_feasible, C, A)
            
            # --- F. 加入种群 ---
            push!(new_population, child_improved)
            
            # 更新全局最好
            if z_child > best_z
                best_z = z_child
                best_x = copy(child_improved)
            end
        end 
        
        # 更新种群
        population = new_population
        
    end # 结束 for 循环
    
    println("  [GA] 结束。最佳 z = ", best_z)
    return best_x, best_z
end

# --------------------------------------------------------------------------- #
#   4. 实验与对比函数
# --------------------------------------------------------------------------- #

"""
    experimentation_EI3()
    运行实验：对比 EI1 (GRASP) 和 EI3 (Advanced GA)
"""
function experimentation_EI3()
    
    println("="^70)
    println("      EI3: BATTLE (GA Improved vs GRASP)")
    println("="^70)
    
    data_directory = "dat"
    # 调用 EI1 的 getfname
    fnames = getfname(data_directory) 
    
    if isempty(fnames)
        return
    end
    
    # 表头
    println("-"^80)
    print(rpad("Instance", 15), " | ")
    print(lpad("EI1(z)", 10), " | ") 
    print(lpad("EI3(GA)", 10), " | ") 
    print(lpad("Imp(%)", 8), " | ")  
    println(lpad("Time(s)", 8))      
    println("-"^80)
    
    for fname in fnames
        fpath = joinpath(data_directory, fname)
        C, A = loadSPP(fpath)
        
        # 1. EI1 基准 (构造 + 局搜)
        x_ei1, z_ei1 = randomized_greedy_construction(C, A, 0.75)
        x_ei1, z_ei1 = local_search(x_ei1, C, A)
        
        # 2. EI3 遗传算法 (2.0秒限制)
        start_ga = time()
        x_ga, z_ga = genetic_algorithm_SPP(C, A, 2.0)
        t_ga = time() - start_ga
        
        # 3. 计算提升
        imp = 0.0
        if z_ei1 > 0
            imp = 100.0 * (z_ga - z_ei1) / z_ei1
        end
        
        # 4. 打印
        s_name = rpad(fname, 15)
        s_ei1  = lpad(string(z_ei1), 10)
        s_ei3  = lpad(string(z_ga), 10)
        s_imp  = lpad(string(round(imp, digits=2)), 8)
        s_t    = lpad(string(round(t_ga, digits=2)), 8)
        
        println(s_name, " | ", s_ei1, " | ", s_ei3, " | ", s_imp, " | ", s_t)
    end
    println("="^80)
end

println("\n'livrableEI3.jl' 已加载。包含贪婪修复和逻辑交叉算子。")

