using JuMP, GLPK
using LinearAlgebra
using Random

# 引用其他写好的文件
# 确保它们都在同一个文件夹里
include("loadSPP.jl")
include("setSPP.jl")
include("getfname.jl")

# 1. 构造解 (GRASP)
# 这里用的是随机贪婪算法
function randomized_greedy_construction(C, A, alpha=0.75)
    m, n = size(A)
    
    # x 是我们的解，一开始全是 0
    x = zeros(Int, n)
    
    # 记录每一行是不是被覆盖了
    row_covered = zeros(Bool, m)
    
    # 一开始所有变量都能选
    candidates = collect(1:n)
    
    while length(candidates) > 0
        
        # 第一步：找出所有还能选的变量
        # (就是它覆盖的行，目前都还没被占用的)
        valid_candidates = Int[]
        for j in candidates
            is_valid = true
            for i = 1:m
                if A[i, j] == 1 && row_covered[i] == true
                    is_valid = false
                    break
                end
            end
            
            if is_valid
                push!(valid_candidates, j)
            end
        end

        # 如果没有变量能选了，就停
        if length(valid_candidates) == 0
            break
        end

        # 第二步：计算阈值 (RCL)
        # 找这些候选人里 C 值最大和最小的
        c_vals = C[valid_candidates]
        c_min = minimum(c_vals)
        c_max = maximum(c_vals)
        
        limit = c_min + alpha * (c_max - c_min)
        
        # 选出那些价值比较高的候选人
        rcl_list = Int[]
        for j in valid_candidates
            if C[j] >= limit
                push!(rcl_list, j)
            end
        end
        
        # 第三步：随机选一个
        if length(rcl_list) == 0
            # 应该不会发生，但为了保险
            j_selected = valid_candidates[1]
        else
            j_selected = rand(rcl_list)
        end

        # 第四步：更新状态
        x[j_selected] = 1
        
        # 标记被这个变量占用的行
        for i = 1:m
            if A[i, j_selected] == 1
                row_covered[i] = true
            end
        end
        
        # 更新候选列表 (把选过的移除掉)
        # 这里用个简单的方法：重建一个列表
        new_candidates = Int[]
        for k in candidates
            if k != j_selected
                push!(new_candidates, k)
            end
        end
        candidates = new_candidates
    end

    # 算一下最后的总价值 z
    z = 0.0
    for j = 1:n
        if x[j] == 1
            z = z + C[j]
        end
    end

    return x, z
end

# 2. 局部搜索 (Local Search)
# 尝试一点点改进我们的解
function local_search(x, C, A)
    m, n = size(A)
    curr_x = copy(x) # 复制一份，别改坏了原来的
    
    # 算出现在的 z 值
    z = 0.0
    for j = 1:n
        if curr_x[j] == 1
            z = z + C[j]
        end
    end
    
    # 只要能找到更好的，就一直循环
    while true
        found_better = false
        
        # --- 方法 1: 尝试加一个变量 (Add) ---
        
        # 打乱顺序来试，这样比较随机
        try_order = shuffle(1:n)
        
        for j in try_order
            # 如果 j 没被选，而且它的价值是正的
            if curr_x[j] == 0 && C[j] > 0
                
                # 检查能不能加 j (是不是会冲突)
                can_add = true
                for i = 1:m
                    if A[i, j] == 1
                        # 看看第 i 行有没有被别人占了
                        for k = 1:n
                            if curr_x[k] == 1 && A[i, k] == 1
                                can_add = false
                                break
                            end
                        end
                    end
                    if can_add == false
                        break
                    end
                end
                
                # 如果可以加，那就加上！
                if can_add
                    curr_x[j] = 1
                    z = z + C[j]
                    found_better = true
                    break # 找到一个就行，重新开始大循环
                end
            end
        end
        
        if found_better
            continue
        end

        # --- 方法 2: 尝试交换 (Swap) ---
        # 找一个在解里的 u，换成一个不在解里的 v
        
        # 先要把“在解里”和“不在解里”的变量分开
        in_sol = Int[]
        out_sol = Int[]
        for j = 1:n
            if curr_x[j] == 1
                push!(in_sol, j)
            else
                push!(out_sol, j)
            end
        end
        
        # 打乱它们
        in_sol = shuffle(in_sol)
        out_sol = shuffle(out_sol)

        for u in in_sol
            for v in out_sol
                # 只有 v 比 u 值钱，才值得换
                if C[v] > C[u]
                    
                    # 检查 v 能不能放进去 (要先忽略 u)
                    can_swap = true
                    
                    for i = 1:m
                        if A[i, v] == 1
                            # 这一行如果有别人占了就不行 (但如果是 u 占的没事，因为 u 要走了)
                            for k = 1:n
                                if curr_x[k] == 1 && k != u
                                    if A[i, k] == 1
                                        can_swap = false
                                        break
                                    end
                                end
                            end
                        end
                        if can_swap == false
                            break
                        end
                    end
                    
                    if can_swap
                        # 执行交换
                        curr_x[u] = 0
                        curr_x[v] = 1
                        z = z - C[u] + C[v]
                        found_better = true
                        break
                    end
                end
            end
            if found_better
                break
            end
        end
        
        # 如果两招都试过了还是没改进，那就结束吧
        if found_better == false
            break
        end
    end

    return curr_x, z
end

# 3. 求最优解 (用 GLPK 求解器)
function solve_optimal(C, A)
    # 建模型
    model = setSPP(C, A)
    set_optimizer(model, GLPK.Optimizer)
    
    # 不显示求解器的那些啰嗦的输出
    set_silent(model) 
    
    t_start = time()
    optimize!(model)
    t_total = time() - t_start
    
    if termination_status(model) == MOI.OPTIMAL
        return objective_value(model), t_total
    else
        return 0.0, t_total
    end
end

# 4. 解决单个文件的函数
function resoudreSPP(fname)
    file_path = "dat/" * fname
    
    # 简单的检查文件存不存在
    if isfile(file_path) == false
        println("找不到文件: " * file_path)
        return
    end

    println("正在处理: " * fname)
    C, A = loadSPP(file_path)

    # 跑我们的算法
    t1 = time()
    x_start, z_start = randomized_greedy_construction(C, A, 0.75)
    x_final, z_final = local_search(x_start, C, A)
    t_heur = time() - t1
    
    println("  我的算法结果: z = ", z_final, "  (用时: ", round(t_heur, digits=3), "秒)")

    # 跑最优解对比一下
    z_opt, t_opt = solve_optimal(C, A)
    println("  最优解结果:   z = ", z_opt, "  (用时: ", round(t_opt, digits=3), "秒)")
    
    if z_opt > 0
        gap = 100 * (z_opt - z_final) / z_opt
        println("  差距 (Gap):   ", round(gap, digits=2), "%")
    end
    println("-"^30)
end

# 5. 批量跑所有测试
function experimentationSPP()
    if isdir("data") == false
        println("错误: 没有 data 文件夹")
        return
    end

    files = getfname("data")
    if length(files) == 0
        println("data 文件夹是空的")
        return
    end

    println("Instance\tHeur\tTime\tOpt\tGap(%)")
    println("-"^50)
    
    for f in files
        C, A = loadSPP("dat/" * f)
        
        t1 = time()
        x0, z0 = randomized_greedy_construction(C, A, 0.75)
        x1, z1 = local_search(x0, C, A)
        t_h = time() - t1
        
        z_opt, t_opt = solve_optimal(C, A)
        
        gap = 0.0
        if z_opt > 0
            gap = 100 * (z_opt - z1) / z_opt
        end
        
        # 打印这一行的结果
        print(f)
        print("\t")
        print(Int(z1))
        print("\t")
        print(round(t_h, digits=3))
        print("\t")
        print(Int(z_opt))
        print("\t")
        println(round(gap, digits=2))
    end
end
