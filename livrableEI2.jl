# =========================================================================== #
#                      livrableEI2.jl (GRASP & Reactive GRASP)
#
#   Exercice d'Implémentation n°2 (EI2)
#   Métaheuristiques - Master Informatique
#
#   Contenu :
#   1. Rappel des fonctions EI1 (Construction, Local Search)
#   2. Métaheuristique GRASP (Classique)
#   3. Métaheuristique Reactive GRASP (Composant additionnel)
#   4. Expérimentation et Comparaison
# =========================================================================== #

println("Chargement des modules...")
using Random
using Statistics # Pour mean()

# On suppose que les fichiers utilitaires sont présents dans le dossier src/
try
    include("loadSPP.jl")
    include("getfname.jl")
catch e
    println("[ERREUR] Fichiers loadSPP.jl ou getfname.jl manquants.")
    println("Assurez-vous d'être dans le dossier 'src'.")
    rethrow(e)
end

# =========================================================================== #
#   PARTIE 1 : RÉUTILISATION DE EI1 (Construction & Recherche Locale)
#   [cite: 26] "Cet exercice ... s'inscrit dans la continuité ... du EI1"
# =========================================================================== #

# --- 1.1 Construction Gloutonne Randomisée (copié/adapté de EI1) ---
function randomized_greedy_construction(C, A, alpha)
    m, n = size(A)
    x = zeros(Int, n)
    solution_indices = Int[]
    covered_rows = zeros(Bool, m)
    available_vars = collect(1:n)

    while !isempty(available_vars)
        # Candidats valides
        candidates = Int[]
        for j in available_vars
            is_valid = true
            for i = 1:m
                if A[i, j] == 1 && covered_rows[i]
                    is_valid = false
                    break
                end
            end
            if is_valid; push!(candidates, j); end
        end

        if isempty(candidates); break; end

        # RCL (Restricted Candidate List)
        C_min = minimum(C[j] for j in candidates)
        C_max = maximum(C[j] for j in candidates)
        threshold = C_min + alpha * (C_max - C_min)

        RCL = [j for j in candidates if C[j] >= threshold]
        if isempty(RCL) 
             # Cas rare de précision flottante
             RCL = [j for j in candidates if C[j] == C_max]
        end

        # Sélection aléatoire
        j_selected = rand(RCL)

        # Mise à jour
        x[j_selected] = 1
        push!(solution_indices, j_selected)
        
        # Marquer les lignes couvertes
        rows_to_cover = [i for i=1:m if A[i, j_selected] == 1]
        for i in rows_to_cover; covered_rows[i] = true; end

        # Mettre à jour available_vars (retirer ceux en conflit)
        # (Version simplifiée beginner : on filtre la liste)
        new_available = Int[]
        for k in available_vars
            if k == j_selected; continue; end
            conflict = false
            for i in rows_to_cover
                if A[i, k] == 1
                    conflict = true
                    break
                end
            end
            if !conflict; push!(new_available, k); end
        end
        available_vars = new_available
    end
    
    z = sum(C[j] for j in solution_indices)
    return x, z
end

# --- 1.2 Recherche Locale (Descente - First Improvement) (copié de EI1) ---
function local_search(x_in, C, A)
    m, n = size(A)
    x = copy(x_in)
    solution_indices = findall(x .== 1)
    z = sum(C[j] for j in solution_indices)
    
    improvement = true
    while improvement
        improvement = false
        
        # 1. Voisinage ADD (Ajout simple)
        candidates_add = shuffle(findall(x .== 0))
        for j in candidates_add
            can_add = true
            for i = 1:m
                if A[i, j] == 1
                    for k in solution_indices
                        if A[i, k] == 1; can_add = false; break; end
                    end
                end
                if !can_add; break; end
            end
            
            if can_add && C[j] > 0
                x[j] = 1
                push!(solution_indices, j)
                z += C[j]
                improvement = true
                break
            end
        end
        if improvement; continue; end

        # 2. Voisinage SWAP (1-1 Echange)
        for idx_out in 1:length(solution_indices)
            j_out = solution_indices[idx_out]
            candidates_in = findall(x .== 0) # On pourrait mélanger ici aussi
            
            for j_in in candidates_in
                if C[j_in] > C[j_out] # Gain positif seulement
                    # Vérifier faisabilité (j_in vs Sol\{j_out})
                    can_swap = true
                    for i = 1:m
                        if A[i, j_in] == 1
                            for k in solution_indices
                                if k != j_out && A[i, k] == 1
                                    can_swap = false; break
                                end
                            end
                        end
                        if !can_swap; break; end
                    end
                    
                    if can_swap
                        x[j_out] = 0
                        x[j_in] = 1
                        solution_indices[idx_out] = j_in
                        z = z - C[j_out] + C[j_in]
                        improvement = true
                        break
                    end
                end
            end
            if improvement; break; end
        end
    end
    return x, z
end

# =========================================================================== #
#   PARTIE 2 : GRASP CLASSIQUE
#    "élaboration d'un solveur reposant sur la métaheuristique GRASP"
# =========================================================================== #

"""
    run_grasp(C, A, max_time, alpha_fixe)

    Exécute un GRASP standard avec un alpha fixe pendant un temps donné.
"""
function run_grasp(C, A, max_time, alpha_fixe)
    start_time = time()
    
    # Meilleure solution trouvée globalement
    best_x = Int[]
    best_z = -Inf
    
    iterations = 0
    
    # Boucle tant qu'il reste du temps [cite: 39]
    while (time() - start_time) < max_time
        iterations += 1
        
        # 1. Construction
        x_curr, z_curr = randomized_greedy_construction(C, A, alpha_fixe)
        
        # 2. Recherche Locale
        x_curr, z_curr = local_search(x_curr, C, A)
        
        # 3. Mise à jour de la meilleure solution
        if z_curr > best_z
            best_z = z_curr
            best_x = copy(x_curr)
        end
    end
    
    return best_x, best_z, iterations
end

# =========================================================================== #
#   PARTIE 3 : REACTIVE GRASP (Composant Additionnel)
#    "(a) ReactiveGRASP"
# =========================================================================== #

"""
    run_reactive_grasp(C, A, max_time)

    Exécute un Reactive GRASP.
    L'alpha n'est pas fixe, il est choisi parmi une liste {0.1, ..., 0.9}.
    Les probabilités de choisir un alpha sont mises à jour périodiquement.
"""
function run_reactive_grasp(C, A, max_time)
    start_time = time()
    
    # Paramètres Reactive
    alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    m_alphas = length(alphas)
    probabilities = fill(1.0/m_alphas, m_alphas) # Au début : équiprobable
    
    # Stats pour chaque alpha : Somme des z, Nombre de fois choisi
    sum_z = zeros(Float64, m_alphas)
    count_chosen = zeros(Int, m_alphas)
    best_z_alpha = zeros(Float64, m_alphas) # Meilleur z trouvé par cet alpha
    
    # Meilleure solution globale
    best_x_global = Int[]
    best_z_global = -Inf
    
    iterations = 0
    block_size = 50 # Mettre à jour les probas toutes les X itérations
    
    # Boucle principale
    while (time() - start_time) < max_time
        iterations += 1
        
        # --- A. Choisir Alpha (Roulette Wheel) ---
        r = rand()
        cumulative_prob = 0.0
        selected_idx = 1
        for i = 1:m_alphas
            cumulative_prob += probabilities[i]
            if r <= cumulative_prob
                selected_idx = i
                break
            end
        end
        alpha_curr = alphas[selected_idx]
        
        # --- B. Exécution GRASP (1 itération) ---
        # 1. Construction
        x_curr, z_curr = randomized_greedy_construction(C, A, alpha_curr)
        # 2. Local Search
        x_curr, z_curr = local_search(x_curr, C, A)
        
        # --- C. Mise à jour des stats ---
        if z_curr > best_z_global
            best_z_global = z_curr
            best_x_global = copy(x_curr)
        end
        
        sum_z[selected_idx] += z_curr
        count_chosen[selected_idx] += 1
        if z_curr > best_z_alpha[selected_idx]
            best_z_alpha[selected_idx] = z_curr
        end
        
        # --- D. Mise à jour des Probabilités (Réactif) ---
        # On fait cela à la fin de chaque "bloc" d'itérations
        if iterations % block_size == 0
            
            # Calculer le score (qualité) de chaque alpha
            # Formule simple : score = z_moyen / z_best_global
            scores = zeros(Float64, m_alphas)
            for i = 1:m_alphas
                if count_chosen[i] > 0
                    avg_z = sum_z[i] / count_chosen[i]
                    # On amplifie les différences (puissance delta, ex: 8)
                    scores[i] = (avg_z / best_z_global)^8 
                else
                    scores[i] = 0.01 # Petite valeur par défaut si jamais choisi
                end
            end
            
            # Normaliser pour obtenir des probabilités (somme = 1)
            total_score = sum(scores)
            if total_score > 0
                probabilities = scores ./ total_score
            else
                probabilities = fill(1.0/m_alphas, m_alphas)
            end
            
            # (Optionnel) Afficher l'alpha favori du moment
            # best_alpha_idx = argmax(probabilities)
            # println("  [Update] Itér $iterations : Alpha favori = ", alphas[best_alpha_idx], " (p=", round(probabilities[best_alpha_idx], digits=2), ")")
        end
        
    end
    
    return best_x_global, best_z_global, iterations, probabilities
end

# =========================================================================== #
#   PARTIE 4 : EXPÉRIMENTATION ET COMPARAISON
#   [cite: 42] "Comparer les résultats obtenus ... avec ... EI1 ... et GRASP"
# =========================================================================== #

function experimentation_EI2()
    println("\n===========================================================")
    println("  DÉMARRAGE DE L'EXPÉRIMENTATION EI2 (Battle GRASP)")
    println("===========================================================")
    
    # Réglages
    time_budget = 30.0 # Secondes par instance [cite: 39] "budget de calcul donné"
    alpha_fixed = 0.75 # Pour le GRASP simple
    
    # Chargement des instances
    fnames = getfname("dat")
    if isempty(fnames)
        println("Aucune instance trouvée dans 'dat/'.")
        return
    end
    
    # Affichage Entête Tableau
    # Nous allons comparer : EI1 vs GRASP(0.75) vs ReactiveGRASP
    println("\n" * "-"^95)
    print(rpad("Instance", 15), " | ")
    print(lpad("EI1(1shot)", 10), " | ")
    print(lpad("GRASP(Fix)", 10), " | ")
    print(lpad("ReactGRASP", 10), " | ")
    print(lpad("BestAlpha", 10), " | ") # Valeur saillante de a [cite: 41]
    println(lpad("NbIter(R)", 9))
    println("-"^95)
    
    for fname in fnames
        fpath = joinpath("dat", fname)
        C, A = loadSPP(fpath)
        
        # 1. EI1 (Construction + LS une seule fois) [cite: 43]
        # On le relance pour être sûr, mais c'est instantané
        x1, z1 = randomized_greedy_construction(C, A, alpha_fixed)
        x1, z1 = local_search(x1, C, A)
        
        # 2. GRASP Classique (alpha = 0.75) [cite: 44]
        # "GRASP sans composant additionnel"
        x_g, z_g, iter_g = run_grasp(C, A, time_budget, alpha_fixed)
        
        # 3. Reactive GRASP (Adaptive) 
        # "Solution opérationnelle"
        x_r, z_r, iter_r, probs = run_reactive_grasp(C, A, time_budget)
        
        # Trouver l'alpha qui a la plus forte probabilité finale
        alphas_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        best_alpha_idx = argmax(probs)
        best_alpha_val = alphas_list[best_alpha_idx]
        
        # Affichage ligne par ligne
        s_name = rpad(fname, 15)
        s_z1   = lpad(string(z1), 10)
        s_zg   = lpad(string(z_g), 10)
        s_zr   = lpad(string(z_r), 10)
        s_alpha= lpad(string(best_alpha_val), 10)
        s_iter = lpad(string(iter_r), 9)
        
        println(s_name, " | ", s_z1, " | ", s_zg, " | ", s_zr, " | ", s_alpha, " | ", s_iter)
        
    end
    println("-"^95)
    println("Budget temps par méthode : ", time_budget, " secondes.")
end

println("\nFichier 'livrableEI2.jl' chargé.")
println("Lancez 'experimentation_EI2()' pour voir la comparaison.")
