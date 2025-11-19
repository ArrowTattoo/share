# =========================================================================== #
#                      livrableEI1.jl (新手友好版)
#
#   这是为 EI1 练习准备的主文件。
#   这个版本被特意编写得更简单、更详细，以便于初学者理解每一步。
#   它优先考虑“易读性”而不是“代码效率”。
# =========================================================================== #

# --------------------------------------------------------------------------- #
#   1. 载入“工具箱” (包) 和 其他 .jl 文件
# --------------------------------------------------------------------------- #

println("正在载入 Julia 的工具箱 (包)...")
using JuMP, GLPK     # 用于数学规划 (解最优解)
using LinearAlgebra  # 这个我们不用了，因为我们会手动计算
using Random         # 用于随机选择 (rand())

println("正在载入你写的其他 .jl 文件...")
try
    # include(...) 就像是把另一个文件的代码复制粘贴到这里
    include("loadSPP.jl")    # 载入 loadSPP 函数
    include("setSPP.jl")     # 载入 setSPP 函数
    include("getfname.jl")   # 载入 getfname 函数
catch e
    println("\n[错误!] 无法找到必要的文件 (如 loadSPP.jl)。")
    println("请确保 livrableEI1.jl 和其他 .jl 文件在同一个 'src' 文件夹中。")
    rethrow(e) # 显示详细的错误信息
end

# --------------------------------------------------------------------------- #
#   2. 启发式构造函数 (随机贪婪算法)
#      [对应 PDF 要求 1]
#
#   目标: 快速“猜”出一个还不错的解
# --------------------------------------------------------------------------- #

"""
    randomized_greedy_construction(C, A, alpha)

    (新手版)
    使用随机贪婪算法 (GRASP) 构造一个解。
    alpha = 0.75 (如同你要求的)
"""
function randomized_greedy_construction(C, A, alpha=0.75)
    
    # m = 约束 (行) 的数量
    # n = 变量 (列) 的数量
    m, n = size(A)
    
    # ----------------
    # 1. 初始化变量
    # ----------------
    
    # 'x' 是我们的解，一开始全是 0 (表示没有选择任何变量)
    x = zeros(Int, n)
    
    # 'solution_indices' 用来存放我们选中的变量的“编号” (比如 [1, 5, 8])
    solution_indices = Int[]
    
    # 'covered_rows' 记录哪些约束行已经被满足了
    # 比如 covered_rows[3] = true 意思是第3行已经有变量了
    covered_rows = zeros(Bool, m)
    
    # 'available_vars' 记录所有“还可以被选择”的变量
    # 一开始，所有变量 (从 1 到 n) 都可以被选择
    available_vars = collect(1:n) # collect() 把 1:n 转换成一个普通数组

    # ----------------
    # 2. 循环构造解
    # ----------------
    
    # 当“可选择”的变量列表不为空时，就继续尝试添加
    while !isempty(available_vars)
        
        # --- a. 找出所有“合法的”候选人 ---
        #    一个候选人 'j' 是合法的，前提是：
        #    它所在的行 (A[i,j]==1) 还没有被满足 (covered_rows[i]==false)
        
        # 'candidates' 存放所有合法候选人的编号
        candidates = Int[]
        
        # 遍历 *每一个* 还在 'available_vars' 列表中的变量
        for j in available_vars
            
            is_valid_candidate = true # 先假设它是合法的
            
            # 检查 'j' 所在的所有行
            for i = 1:m
                # 如果 1) 'j' 在第 i 行 (A[i,j]==1) 
                # 并且 2) 第 i 行已经被占用了 (covered_rows[i] == true)
                if A[i, j] == 1 && covered_rows[i] == true
                    
                    # 那么 'j' 就不是一个合法的候选人
                    is_valid_candidate = false
                    break # 跳出这个内层循环 (不用再检查 'j' 的其他行了)
                end
            end
            
            # 如果 'j' 通过了所有行的检查
            if is_valid_candidate
                # 把它加入到“候选人”列表
                push!(candidates, j)
            end
        end # 结束 'j' 的遍历

        # --- b. 检查是否卡住了 ---
        if isempty(candidates)
            # 如果 'candidates' 列表是空的，说明没有变量能再加入了
            # 我们就完成了构造
            break # 跳出 'while' 循环
        end

        # --- c. 建立“受限候选列表 (RCL)” ---
        #    这是 GRASP 算法的核心：不是总选最好的，而是选“比较好”的
        
        # 1. 找到候选人中的“最好” (C_max) 和“最差” (C_min) 的 C 值
        C_min = Inf    # 设为正无穷大
        C_max = -Inf   # 设为负无穷大
        
        for j in candidates
            if C[j] > C_max
                C_max = C[j]
            end
            if C[j] < C_min
                C_min = C[j]
            end
        end
        
        # 2. 用你指定的 alpha=0.75 计算“阈值”
        threshold = C_min + alpha * (C_max - C_min)

        # 3. 建立 RCL 列表
        RCL = Int[] # RCL = Restricted Candidate List
        for j in candidates
            # 如果这个候选人的 C 值“足够好” (大于等于阈值)
            if C[j] >= threshold
                # 就把它加入 RCL
                push!(RCL, j)
            end
        end
        
        # 4. (安全检查) 如果因为计算问题 RCL 是空的,
        #    我们就只把最好的那个加进去
        if isempty(RCL)
            for j in candidates
                if C[j] == C_max
                    push!(RCL, j)
                end
            end
        end


        # --- d. 从 RCL 中“随机”选一个 ---
        j_selected = rand(RCL) # rand() 会从列表里随机挑一个

        # --- e. 更新我们的解 ---
        
        # 1. 把选中的 'j_selected' 加入解
        x[j_selected] = 1
        push!(solution_indices, j_selected)
        
        # 2. 标记 'j_selected' 占用的“行”
        rows_covered_by_j = Int[] # 记录它占了哪些行
        for i = 1:m
            if A[i, j_selected] == 1
                covered_rows[i] = true
                push!(rows_covered_by_j, i) # 记住第 i 行
            end
        end
        
        # 3. 更新 'available_vars' 列表
        #    (这是最复杂的一步)
        #    我们必须移除 'j_selected' 
        #    以及 *所有* 与 'j_selected' 冲突的变量 (即碰了同一行的变量)
        
        vars_to_remove = Int[] # 建立一个“待移除”列表
        push!(vars_to_remove, j_selected) # 1. 移除 j_selected 自己
        
        # 2. 找出所有冲突的变量
        for i in rows_covered_by_j # 遍历 j 占用的行
            for k in available_vars   # 遍历所有还可用的变量
                if A[i, k] == 1 # 如果 k 也在第 i 行
                    # 那么 k 也必须被移除
                    push!(vars_to_remove, k)
                end
            end
        end
        
        # 3. (新手的方式) 创建一个新的 available_vars 列表
        new_available_vars = Int[]
        for j in available_vars # 遍历旧的列表
            
            # 检查 'j' 是否在 'vars_to_remove' 列表里
            should_remove = false
            for r in vars_to_remove
                if j == r
                    should_remove = true
                    break
                end
            end
            
            # 如果 'j' *不在* 待移除列表里
            if !should_remove
                # 就把它加入到 *新* 列表
                push!(new_available_vars, j)
            end
        end
        
        # 4. 用新列表替换旧列表
        available_vars = new_available_vars
        
    end # 结束 'while' 循环

    # ----------------
    # 3. 计算最终的目标函数值 z
    # ----------------
    z = 0.0 # 用 0.0 来确保结果是浮点数
    
    # 遍历我们选中的所有变量 (j_sol)
    for j_sol in solution_indices
        # 把它的 C 值加到总分 z
        z = z + C[j_sol]
    end
    
    # 另一种方法 (遍历 x 向量):
    # for j = 1:n
    #     if x[j] == 1
    #         z = z + C[j]
    #     end
    # end
    
    # 返回构造好的解 'x' 和它的值 'z'
    return x, z
end


# --------------------------------------------------------------------------- #
#   3. 局部搜索函数 (Descente)
#      [对应 PDF 要求 2]
#
#   目标: 拿到构造的解，尝试“微调”它，让它变得更好
# --------------------------------------------------------------------------- #

"""
    local_search(x_in, C, A)

    (新手版)
    尝试改进一个给定的解 'x_in'。
    我们用两种方法 (邻域) 尝试改进：
    N1: Add (添加一个变量)
    N2: Swap (交换一个变量)
"""
function local_search(x_in, C, A)
    m, n = size(A)
    
    # 1. 复制解 (我们不想修改原始的 'x_in')
    x = copy(x_in)
    
    # 2. 找到这个解的当前 z 值
    z = 0.0
    solution_indices = Int[] # 同时，找到所有 x[j]==1 的 j
    
    for j = 1:n
        if x[j] == 1
            z = z + C[j]
            push!(solution_indices, j) # 把 j 存起来
        end
    end
    
    # 3. 循环改进
    improvement_found = true # 这是一个“开关”
    
    # 只要上一轮找到了改进 (开关=true), 就再试一轮
    while improvement_found == true
        
        # 每一轮开始时，先把“开关”关掉
        improvement_found = false 
        
        # --- 邻域 N1 : ADD (尝试添加一个变量) ---
        
        # 遍历所有变量
        for j_in = 1:n
            
            # 如果 j_in 还不在解里 (x[j_in] == 0)
            if x[j_in] == 0
                
                # 检查 'j_in' 能不能被安全地加入
                can_add = true # 先假设可以
                
                for i = 1:m # 遍历所有约束行
                    # 如果 j_in 在第 i 行 (A[i, j_in] == 1)
                    if A[i, j_in] == 1
                        # 我们必须检查第 i 行是否已被占用
                        for k_sol in solution_indices # 遍历已在解中的 k
                            if A[i, k_sol] == 1
                                # 冲突！k_sol 也在第 i 行
                                can_add = false
                                break # 停止检查 k_sol
                            end
                        end
                    end
                    if !can_add; break; end # 停止检查 i
                end
                
                # 如果 1) 可以添加 并且 2) 这个变量是“有价值的” (C[j_in] > 0)
                if can_add && C[j_in] > 0
                    # 太好了！找到了一个改进！
                    x[j_in] = 1         # 1. 更新 x
                    push!(solution_indices, j_in) # 2. 更新解列表
                    z = z + C[j_in]     # 3. 更新 z 值
                    
                    improvement_found = true # 4. 打开“开关”
                    
                    # println("LS (Add): 找到了! +", j_in) # (用于调试)
                    
                    # 采用“第一个改进”策略：一旦找到，马上停止搜索
                    break 
                end
            end
        end # 结束 N1 的 j_in 循环
        
        # 如果 N1 找到了改进 (开关=true),
        # 我们就 'continue' (跳过 N2, 直接开始下一轮 while 循环)
        if improvement_found == true
            continue 
        end
        
        # --- 邻域 N2 : SWAP (尝试 1-换-1) ---
        #    用一个“局外”变量 (k_in) 替换一个“局内”变量 (j_out)
        
        # 遍历所有“局内”变量
        for j_out in solution_indices
            
            # 遍历所有“局外”变量
            for k_in = 1:n
                if x[k_in] == 0 # k_in 必须是“局外”的
                    
                    # 1. 检查这个交换是否“划算”
                    gain = C[k_in] - C[j_out]
                    if gain > 0 # 必须是正收益
                        
                        # 2. 检查这个交换是否“合法”
                        #    (即 k_in 是否会和 j_out *以外* 的变量冲突)
                        can_swap = true # 先假设合法
                        
                        for i = 1:m # 遍历所有约束行
                            # 如果 k_in 在第 i 行
                            if A[i, k_in] == 1
                                # 检查 k_in 是否和“其他”解冲突
                                for l_sol in solution_indices
                                    
                                    # l_sol == j_out ? 我们不在乎
                                    # (我们正要把 j_out 换出去)
                                    if l_sol != j_out 
                                        if A[i, l_sol] == 1
                                            # 冲突！k_in 和 l_sol 都在第 i 行
                                            can_swap = false
                                            break
                                        end
                                    end
                                end # 结束 l_sol 循环
                            end
                            if !can_swap; break; end # 结束 i 循环
                        end
                        
                        # 如果 1) 划算 并且 2) 合法
                        if can_swap
                            # 太好了！又找到了一个改进！
                            x[j_out] = 0        # 1. 移除 j_out
                            x[k_in] = 1         # 2. 加入 k_in
                            z = z + gain        # 3. 更新 z 值
                            
                            # 4. 更新解列表 (这是最棘手的一步)
                            #    我们要把 j_out 替换成 k_in
                            
                            # 先找到 j_out 在列表中的“位置”
                            idx_to_remove = -1
                            for i = 1:length(solution_indices)
                                if solution_indices[i] == j_out
                                    idx_to_remove = i
                                    break
                                end
                            end
                            
                            # 在那个位置换上 k_in
                            if idx_to_remove != -1
                                solution_indices[idx_to_remove] = k_in
                            end
                            
                            improvement_found = true # 5. 打开“开关”
                            
                            # println("LS (Swap): 找到了! ", j_out, "->", k_in) # (用于调试)
                            
                            # 采用“第一个改进”策略
                            break # 停止 k_in 循环
                        end
                    end
                end
            end # 结束 k_in 循环
            
            if improvement_found == true
                break # 停止 j_out 循环
            end
        end # 结束 N2 的 j_out 循环
        
        # 如果 N1 和 N2 都没找到改进 (开关=false),
        # 那么 while 循环将在下次检查时停止
        
    end # 结束 'while improvement_found' 循环
    
    return x, z
end


# --------------------------------------------------------------------------- #
#   4. 精确求解函数 (用 JuMP + GLPK)
#      [对应 PDF 要求 4]
#
#   目标: 调用专业求解器，计算出“真正的”最优解
# --------------------------------------------------------------------------- #

"""
    solve_optimal(C, A)

    使用 JuMP 和 GLPK 求解器计算 SPP 的精确最优解。
"""
function solve_optimal(C, A)
    
    # 1. 调用 setSPP.jl 里的函数来建立数学模型
    #    (spp_model 是一个 JuMP 模型对象)
    spp_model = setSPP(C, A)
    
    # 2. 告诉 JuMP 我们想用哪个求解器
    #    (我们用 GLPK，一个免费的开源求解器)
    set_optimizer(spp_model, GLPK.Optimizer)
    
    # 3. 开始计时
    start_time = time()
    
    # 4. 运行求解器！
    println("  (GLPK 正在计算最优解...)")
    optimize!(spp_model)
    
    # 5. 停止计时
    t_optimal_seconds = time() - start_time
    
    # 6. 检查求解器是否成功
    z_optimal = -1.0 # 默认值
    status = termination_status(spp_model)
    
    if status == MOI.OPTIMAL
        # 成功！获取最优解的值
        z_optimal = objective_value(spp_model)
    else
        # 失败！
        println("[警告] GLPK 求解器未能找到最优解。状态: ", status)
    end
    
    return z_optimal, t_optimal_seconds
end


# --------------------------------------------------------------------------- #
#   5. 函数 `resoudreSPP(fname)`
#      [PDF (第2页) 要求 2]
#
#   目标: 解决 *一个* 实例，并打印所有结果
# --------------------------------------------------------------------------- #

"""
    resoudreSPP(fname)

    (主函数1)
    加载、解决 (启发式 + 精确解) 并打印单个实例的结果。
    'fname' 应该是文件名, 比如 "didactic.dat"
"""
function resoudreSPP(fname)
    
    println("="^50) # 打印 50 个 = 作为分隔线
    println("  开始处理单个实例 : ", fname)
    println("="^50)
    
    # --- 1. 加载数据 ---
    # PDF 要求文件在 'dat' 文件夹下
    fpath = joinpath("dat", fname) 
    # joinpath 会聪明地组合路径 (比如 "dat/didactic.dat")
    
    # 检查文件是否存在
    if !isfile(fpath)
        println("[错误!] 找不到文件: ", fpath)
        println("请确保你的实例文件在 'dat' 子文件夹中。")
        return # 停止函数
    end
    
    C, A = loadSPP(fpath)
    m, n = size(A)
    println("实例加载成功: ", m, " 个约束 (行), ", n, " 个变量 (列)")
    
    # --- 2. 运行我们的启发式算法 ---
    println("\n--- 步骤 1: 运行启发式算法 ---")
    
    # 开始计时
    start_heuristic_time = time()
    
    # a. 构造
    x_const, z_const = randomized_greedy_construction(C, A, 0.75)
    
    # b. 局部搜索 (用构造的解 x_const 作为起点)
    x_local, z_local = local_search(x_const, C, A)
    
    # 停止计时
    t_heuristic_total = time() - start_heuristic_time
    
    println("  构造 (α=0.75) 结果 : z = ", z_const)
    println("  局部搜索结果          : z = ", z_local)
    println("  启发式算法总耗时       : ", round(t_heuristic_total, digits=4), " 秒")
    
    
    # --- 3. 运行精确求解器 ---
    println("\n--- 步骤 2: 运行精确求解器 (GLPK) ---")
    
    z_opt, t_opt = solve_optimal(C, A) # [对应 PDF 要求 4]
    
    println("  GLPK 找到的最优解     : z = ", z_opt)
    println("  精确求解总耗时         : ", round(t_opt, digits=4), " 秒")
    
    # --- 4. 比较和讨论 ---
    println("\n--- 步骤 3: 结果对比 ---")
    
    if z_opt > 0
        # 计算“差距”(Gap)
        # Gap = (最优解 - 我们的解) / 最优解
        gap_percent = 100.0 * (z_opt - z_local) / z_opt
        
        println("  我们的解 (", z_local, ") 达到了最优解 (", z_opt, ") 的 ", round(100.0 * z_local / z_opt, digits=2), "%")
        println("  -> 差距 (Gap) : ", round(gap_percent, digits=2), "%")
    end
    
    println("  速度对比 (启发式 vs 精确) : ", round(t_heuristic_total, digits=4), "s  vs  ", round(t_opt, digits=4), "s")
    println("-"^50, "\n")
end


# --------------------------------------------------------------------------- #
#   6. 函数 `experimentationSPP()`
#      [PDF (第2页) 要求 3]
#
#   目标: 解决 'dat' 文件夹下的 *所有* 实例，并打印总结表格
# --------------------------------------------------------------------------- #

"""
    experimentationSPP()

    (主函数2)
    在 "dat/" 目录下的所有实例上运行实验。
"""
function experimentationSPP()
    
    println("="^70)
    println("      开始批量实验 (EI1)")
    println("="^70)
    
    # PDF 要求数据目录为 'dat'
    data_directory = "dat" 
    
    # --- 1. 获取所有实例的文件名 ---
    println("正在 '", data_directory, "' 文件夹中搜索实例...")
    fnames = [] # 准备一个空列表
    try
        # 调用 getfname.jl 里的函数
        fnames = getfname(data_directory)
        
        if isempty(fnames)
             println("[警告] 在 '", data_directory, "' 文件夹中没有找到任何实例文件。")
             println("请将你的 .dat 文件放入 'dat' 子文件夹。")
             return # 停止函数
        end
        
        println("找到了 ", length(fnames), " 个实例: ", fnames)
        
        # [PDF 要求 "au moins 10 instances" (至少10个)]
        if length(fnames) < 10
            println("[提示] 你找到了 ", length(fnames), " 个实例。PDF 要求至少 10 个。")
        end

    catch e
        println("\n[错误!] 无法读取 '", data_directory, "' 文件夹。")
        println("请确保你创建了 '", data_directory, "' 子文件夹。")
        return
    end

    # --- 2. 准备空列表来存储所有结果 ---
    results_name = String[]      # 存文件名
    results_m = Int[]            # 存 m (行数)
    results_n = Int[]            # 存 n (列数)
    results_z_local = Float64[]  # 存我们的解
    results_t_heuristic = Float64[] # 存我们的时间
    results_z_optimal = Float64[]   # 存最优解
    results_t_optimal = Float64[]   # 存最优解的时间
    
    println("\n--- 开始循环测试所有实例 ---")
    
    # --- 3. 循环遍历每个文件名 ---
    for (idx, fname) in enumerate(fnames) # enumerate 给我们 索引(idx) 和 值(fname)
        
        print("  (", idx, "/", length(fnames), ") 正在处理: ", fname, "...")
        
        # 加载
        fpath = joinpath(data_directory, fname)
        C, A = loadSPP(fpath)
        m, n = size(A)
        
        # 运行启发式
        start_h = time()
        x_const, z_const = randomized_greedy_construction(C, A, 0.75)
        x_local, z_local = local_search(x_const, C, A)
        t_h = time() - start_h
        
        # 运行精确解
        z_opt, t_opt = solve_optimal(C, A)
        
        # 保存结果
        push!(results_name, fname)
        push!(results_m, m)
        push!(results_n, n)
        push!(results_z_local, z_local)
        push!(results_t_heuristic, t_h)
        push!(results_z_optimal, z_opt)
        push!(results_t_optimal, t_opt)
        
        println(" 完成.")
        
    end # 结束文件循环
    
    println("--- 所有实例处理完毕 ---")
    
    # --- 4. 打印最终的总结表格 ---
    println("\n" * "="^70)
    println("                    实验结果总结表格")
    println("="^70)
    
    # 打印表头
    # 我们用 lpad (left-pad, 向左填充空格) 和 rpad (right-pad) 来对齐
    print(rpad("实例", 15), " | ")
    print(lpad("m", 5), " | ")
    print(lpad("n", 5), " | ")
    print(lpad("z_Local", 10), " | ")
    print(lpad("t_Heur(s)", 9), " | ")
    print(lpad("z_Optimal", 10), " | ")
    print(lpad("t_Opt(s)", 9), " | ")
    println(lpad("Gap(%)", 7)) # println 会换行
    
    println("-"^79) # 打印分隔线
    
    # 循环打印每一行的数据
    for i = 1:length(fnames)
        
        # 计算 Gap
        gap = 0.0
        if results_z_optimal[i] > 0
            gap = 100.0 * (results_z_optimal[i] - results_z_local[i]) / results_z_optimal[i]
        end
        
        # 把数字转换成格式化的字符串
        s_name  = rpad(results_name[i], 15)
        s_m     = lpad(string(results_m[i]), 5)
        s_n     = lpad(string(results_n[i]), 5)
        s_z_loc = lpad(string(round(results_z_local[i], digits=1)), 10)
        s_t_heu = lpad(string(round(results_t_heuristic[i], digits=4)), 9)
        s_z_opt = lpad(string(round(results_z_optimal[i], digits=1)), 10)
        s_t_opt = lpad(string(round(results_t_optimal[i], digits=4)), 9)
        s_gap   = lpad(string(round(gap, digits=2)), 7)
        
        # 打印
        println(s_name, " | ", s_m, " | ", s_n, " | ", s_z_loc, " | ", s_t_heu, " | ", s_z_opt, " | ", s_t_opt, " | ", s_gap)
    end
    
    println("="^70)
    println("实验结束。")
end

# --------------------------------------------------------------------------- #
#   文件末尾
# --------------------------------------------------------------------------- #

println("\n'livrableEI1.jl' (新手版) 已经加载完毕。")
println("你可以开始运行了 (在 Julia REPL 中):")
println("  -> resoudreSPP(\"didactic.dat\")  (注意：'dat' 文件夹中必须有这个文件)")
println("  -> experimentationSPP()")
