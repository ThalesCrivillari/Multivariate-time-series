# ME731 - Projeto P1 - RA 236312 - PNAR aplicado a dengue (RMC, 2015-2024)
# Gera as figuras em figuras/ e imprime no console os numeros usados no relatorio.

#pacotes
install.packages(c("PNAR", "sf", "spdep", "igraph", "surveillance", "ggplot2",
                   "patchwork", "dplyr", "tidyr", "readr", "tseries"))
library(PNAR)
library(sf)
library(spdep)
library(igraph)
library(surveillance)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(tseries)

#configuracoes
sf_use_s2(FALSE)
fig = function(g, nome, w, h) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(file.path("figuras", nome), g, device = cairo_pdf, width = w, height = h)
}
tab = function(titulo, df, ...) {
  cat("\n##", titulo, "\n")
  print(as.data.frame(df), row.names = FALSE)
}


#leitura dos dados: tudo lido do GitHub (copia dos dados originais)
#fontes originais: malha municipal do IBGE (2025), sedes municipais (geobr/IPEA, 2010) e casos semanais do InfoDengue (Fiocruz/FGV)
repo = "https://raw.githubusercontent.com/ThalesCrivillari/Multivariate-time-series/main/dados/"

#malha dos municipios (o shapefile tem 4 arquivos, baixados numa pasta temporaria)
pasta = tempdir()
for (ext in c("shp", "shx", "dbf", "prj")) {
  arquivo = paste0("SP_Municipios_2025.", ext)
  download.file(paste0(repo, arquivo), file.path(pasta, arquivo), mode = "wb", quiet = TRUE)
}
malha_sp = st_read(file.path(pasta, "SP_Municipios_2025.shp"), quiet = TRUE) %>%
  transmute(cod = as.numeric(CD_MUN), nome = NM_MUN, rgint = NM_RGINT)
malha_sp = st_make_valid(st_transform(malha_sp, 4326))
malha_plot = st_transform(st_simplify(st_transform(malha_sp, 5880), preserveTopology = TRUE, dTolerance = 400), 4326)

#coordenadas das sedes
sedes = read_csv(paste0(repo, "sedes_sp.csv"), show_col_types = FALSE)
xy_sp = as.matrix(sedes[match(malha_sp$cod, sedes$geocode), c("lon", "lat")])

#municipios: RMC (20) e regiao ampliada (87)
cod_rmc = c(3501608, 3503802, 3509502, 3512803, 3515152, 3519055, 3519071, 3520509,
            3523404, 3524709, 3531803, 3532009, 3533403, 3536505, 3537107, 3545803,
            3548005, 3552403, 3556206, 3556701)
cod_grande = malha_sp$cod[malha_sp$rgint %in% "Campinas"]

#casos semanais de dengue (um arquivo por municipio, colunas geocode, SE, casos, pop)
dengue = bind_rows(lapply(cod_grande, function(g) {
  read_csv(paste0(repo, "infodengue/", g, ".csv"), show_col_types = FALSE)
})) %>%
  filter(SE >= 201501, SE <= 202452) %>%
  distinct(geocode, SE, .keep_all = TRUE)

#funcoes: distancias, redes e matrizes W
distancia = function(lon, lat) {
  f0 = mean(lat) * pi / 180
  as.matrix(dist(cbind(111.195 * cos(f0) * lon, 111.195 * lat)))
}
normaliza = function(A) {
  W = A / rowSums(A)
  W[is.na(W)] = 0
  W
}

prepara = function(cods) {
  cods = sort(cods)
  ii = match(cods, malha_sp$cod)
  Yw = dengue %>% filter(geocode %in% cods) %>% select(SE, geocode, casos) %>%
    pivot_wider(names_from = geocode, values_from = casos) %>% arrange(SE)
  SE = Yw$SE
  Y = as.matrix(Yw[, as.character(cods)])
  Y[is.na(Y)] = 0
  colnames(Y) = malha_sp$nome[ii]
  pop = dengue %>% filter(geocode %in% cods) %>% group_by(geocode) %>%
    summarise(pop = median(pop, na.rm = TRUE)) %>% arrange(geocode) %>% pull(pop)
  stopifnot(!anyNA(pop))
  xy = xy_sp[ii, ]
  D = distancia(xy[, 1], xy[, 2])
  diag(D) = Inf
  malha = malha_sp[ii, ]
  A_cont = nb2mat(poly2nb(malha, queen = TRUE), style = "B", zero.policy = TRUE)
  A_knn = t(apply(D, 1, function(d) as.numeric(rank(d, ties.method = "first") <= 4)))
  dstar = max(apply(D, 1, min))
  A_dist = ifelse(D <= dstar, 1 / D, 0)
  A_grav = sweep(1 / D^2, 2, pop, "*")
  diag(A_cont) = diag(A_knn) = diag(A_dist) = diag(A_grav) = 0
  redes = list("Contiguidade" = A_cont, "k-vizinhos (k=4)" = A_knn,
                "Distância inversa" = A_dist, "Gravitacional" = A_grav)
  treino = SE < 202401
  list(cods = cods, Y = Y, SE = SE, pop = pop, share = pop / sum(pop), xy = xy,
       malha = malha, malha_plot = malha_plot[ii, ], redes = redes, Ws = lapply(redes, normaliza), dstar = dstar,
       treino = treino, Tn = sum(treino), Ytr = Y[treino, ], N = length(cods), TT = nrow(Y),
       datas = as.Date(paste0(substr(SE, 1, 4), "-01-01")) + (SE %% 100 - 1) * 7)
}
tabela_redes = function(r) bind_rows(lapply(names(r$redes), function(n) {
  B = r$redes[[n]] > 0
  g = graph_from_adjacency_matrix(B * 1, mode = "directed")
  setNames(data.frame(n, sum(B), mean(B[row(B) != col(B)]), mean(rowSums(B)),
                      min(rowSums(B)), components(g, mode = "weak")$no),
           c("Rede", "Ligações", "Densidade", "Grau médio",
             "Grau mínimo", "Componentes"))
}))
#funcoes: ajuste do PNAR (ordem p, rede, AIC/BIC/QIC)
ajusta_grade = function(r, ps = 1:8, p_log = 2) {
  pmax_ = max(ps)
  Zl = list("Linear" = matrix(r$share, ncol = 1), "Log-linear" = matrix(log(r$share), ncol = 1))
  grade = expand.grid(rede = names(r$Ws), tipo = c("Linear", "Log-linear"), p = ps,
                       stringsAsFactors = FALSE)
  grade$AIC = grade$BIC = grade$QIC = grade$soma = NA_real_
  for (rede in names(r$Ws)) {
    ini = NULL
    for (tipo in c("Linear", "Log-linear")) for (p in ps) {
      Y = r$Ytr[(pmax_ - p + 1):nrow(r$Ytr), ]
      f = if (tipo == "Linear") {
        lin_estimnarpq(Y, r$Ws[[rede]], p, Z = Zl[[tipo]], uncons = TRUE, init = ini, maxeval = 2000)
      } else {
        log_lin_estimnarpq(Y, r$Ws[[rede]], p, Z = Zl[[tipo]], uncons = TRUE, maxeval = 1000)
      }
      b = f$coefs[, 1]
      if (tipo == "Linear") ini = c(b[1], b[2:(p + 1)], 1e-4, b[(p + 2):(2 * p + 1)], 1e-4, b[length(b)])
      j = which(grade$rede == rede & grade$tipo == tipo & grade$p == p)
      grade[j, c("AIC", "BIC", "QIC")] = f$ic[c("AIC", "BIC", "QIC")]
      grade$soma[j] = sum(abs(b[2:(2 * p + 1)]))
    }
  }
  melhor = bind_rows(
    grade %>% filter(tipo == "Linear") %>% slice_min(QIC, n = 1),
    grade %>% filter(tipo == "Log-linear", p == p_log) %>% slice_min(QIC, n = 1))
  best = list()
  for (i in seq_len(nrow(melhor))) {
    tipo = melhor$tipo[i]
    p = melhor$p[i]
    W = r$Ws[[melhor$rede[i]]]
    f = if (tipo == "Linear") {
      lin_estimnarpq(r$Ytr, W, p, Z = Zl[[tipo]], uncons = TRUE, maxeval = 2000)
    } else {
      log_lin_estimnarpq(r$Ytr, W, p, Z = Zl[[tipo]], uncons = TRUE, maxeval = 1000)
    }
    best[[tipo]] = list(fit = f, W = W, p = p, rede = melhor$rede[i], Z = Zl[[tipo]])
  }
  list(grade = grade, melhor = melhor, best = best)
}

#funcoes: media estimada, tabela de coeficientes, RPS e previsao
lambda_pnar = function(obj, tipo, Yall, t_idx) {
  b = obj$fit$coefs[, 1]
  p = obj$p
  W = obj$W
  f = if (tipo == "Linear") identity else log1p
  t(sapply(t_idx, function(t) {
    eta = b[1] + b[2 * p + 2] * obj$Z[, 1]
    for (h in 1:p) {
      x = f(Yall[t - h, ])
      eta = eta + b[1 + h] * as.vector(W %*% x) + b[1 + p + h] * x
    }
    if (tipo == "Linear") eta else exp(eta)
  }))
}

tabela_coef = function(best) bind_rows(lapply(names(best), function(m) {
  cf = best[[m]]$fit$coefs
  par_tex = rownames(cf)
  setNames(data.frame(m, par_tex, cf[, 1], cf[, 2], cf[, 3], cf[, 4]),
           c("Modelo", "Parâmetro", "Estimativa", "EP", "z", "p-valor"))
}))

rps_pois = function(y, mu) {
  mapply(function(yy, m) {
    if (!is.finite(m) || m > 1e5) return(abs(yy - min(m, 1e12)))
    k = 0:max(yy, qpois(1 - 1e-10, m))
    sum((ppois(k, m) - (yy <= k))^2)
  }, c(y), c(mu))
}

preve = function(r, best, extra = list()) {
  idx = which(!r$treino)
  Y = r$Y
  Yte = Y[idx, ]
  prev = list()
  for (m in names(best)) prev[[paste("PNAR", m)]] = lambda_pnar(best[[m]], m, Y, idx)
  Xlag = log1p(Y[-r$TT, ])
  prev[["VAR Poisson irrestrito"]] = suppressWarnings(sapply(1:r$N, function(i) {
    df = data.frame(y = Y[-1, i], Xlag)
    g = glm(y ~ ., family = poisson, data = df[1:(r$Tn - 1), ])
    predict(g, newdata = df[idx - 1, ], type = "response")
  }))
  prev[["Ingênuo"]] = pmax(Y[idx - 1, ], 0.5)
  for (m in names(extra)) prev[[m]] = extra[[m]]
  tab = bind_rows(lapply(names(prev), function(m) {
    mu = prev[[m]]
    rps = mean(rps_pois(Yte, mu))
    data.frame(Modelo = m, MAE = mean(abs(Yte - mu)), RMSE = sqrt(mean((Yte - mu)^2)), RPS = rps)
  }))
  list(prev = prev, tab = tab[order(tab$RPS), ], idx = idx)
}

#funcoes: extensao sazonal (termos de Fourier)
pnar_sazonal = function(r, W, p, K) {
  ly = log1p(r$Y)
  wly = ly %*% t(W)
  sem = r$SE %% 100
  monta = function(tt) {
    X = cbind(1, do.call(cbind, lapply(1:p, function(h) c(t(wly[tt - h, , drop = FALSE])))),
               do.call(cbind, lapply(1:p, function(h) c(t(ly[tt - h, , drop = FALSE])))),
               rep(log(r$share), length(tt)))
    if (K > 0) for (k in 1:K)
      X = cbind(X, rep(sin(2 * pi * k * sem[tt] / 52), each = r$N),
                    rep(cos(2 * pi * k * sem[tt] / 52), each = r$N))
    colnames(X) = c("beta0", paste0("beta1", 1:p), paste0("beta2", 1:p), "delta1",
                     if (K > 0) paste0(rep(c("sen", "cos"), K), rep(1:K, each = 2)))
    X
  }
  tt = (p + 1):r$Tn
  X = monta(tt)
  y = c(t(r$Y[tt, ]))
  fit = glm.fit(X, y, family = poisson())
  b = fit$coefficients
  lam = as.vector(exp(X %*% b))
  H = crossprod(X * lam, X)
  S = rowsum(X * (y - lam), rep(tt, each = r$N))
  B = crossprod(S)
  V = solve(H) %*% B %*% solve(H)
  ep = sqrt(diag(V))
  list(coefs = data.frame(Estimativa = b, EP = ep, z = b / ep, p = 2 * pnorm(-abs(b / ep)),
                          row.names = colnames(X)),
       QIC = -2 * sum(y * log(lam) - lam) + 2 * sum(diag(B %*% solve(H))),
       soma = sum(abs(b[2:(2 * p + 1)])), p = p, K = K, b = b, monta = monta)
}
prev_sazonal = function(obj, r, idx) {
  matrix(exp(obj$monta(idx) %*% obj$b), nrow = length(idx), byrow = TRUE)
}

#validacao: exemplo do pacote (influenza)
data(fluBYBW)
fit_flu = lin_estimnarpq(observed(fluBYBW), (neighbourhood(fluBYBW) == 1) * 1, p = 2,
                          Z = matrix(population(fluBYBW)[1, ], ncol = 1), maxeval = 1000)
tab("Influenza (pacote)", data.frame(par = rownames(fit_flu$coefs), fit_flu$coefs[, 1:2]), 4)

#redes e identificacao (RMC)
r = prepara(cod_rmc)
Y = r$Y
Ytr = r$Ytr
N = r$N
Tn = r$Tn
pop = r$pop
datas = r$datas
cod = r$cods
tab("Redes (RMC)", tabela_redes(r), 2)
cat(sprintf("d* = %.1f km\n", r$dstar))

#imagens: redes
arestas = bind_rows(lapply(names(r$Ws), function(n) {
  W = r$Ws[[n]]
  id = which(W > 0, arr.ind = TRUE)
  data.frame(rede = factor(n, names(r$Ws)), x = r$xy[id[, 1], 1], y = r$xy[id[, 1], 2],
             xend = r$xy[id[, 2], 1], yend = r$xy[id[, 2], 2], w = W[id])
}))
g_redes = ggplot() + theme_bw(base_size = 9) +
  geom_sf(data = r$malha_plot, fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_segment(data = arestas, aes(x, y, xend = xend, yend = yend, linewidth = w, alpha = w), colour = "#1f5a99") +
  geom_point(data = data.frame(lon = r$xy[, 1], lat = r$xy[, 2], pop = pop), aes(lon, lat, size = pop), colour = "#b2182b") +
  scale_alpha(range = c(0.35, 0.9), guide = "none") + scale_linewidth(range = c(0.35, 1.6), guide = "none") +
  scale_size(range = c(0.6, 3.5), guide = "none") + facet_wrap(~rede, ncol = 2) + labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
fig(g_redes, "redes.pdf", 6.5, 6.2)

#imagens: series e correlacao
long = data.frame(data = rep(datas, N), mun = rep(colnames(Y), each = r$TT), y = c(Y))
ordem = colnames(Y)[order(pop)]
g_calor = ggplot(long, aes(data, factor(mun, levels = ordem), fill = log1p(y))) + theme_bw(base_size = 9) + geom_tile() +
  scale_fill_viridis_c(name = "log(1+Y)") + scale_x_date(expand = c(0, 0)) + labs(x = NULL, y = NULL)
g_total = ggplot(data.frame(data = datas, total = rowSums(Y)), aes(data, total)) + theme_bw(base_size = 9) + geom_line(linewidth = 0.3) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
  geom_vline(xintercept = as.Date("2024-01-01"), linetype = 2, colour = "red") + labs(x = NULL, y = "Casos na RMC (escala log)")
g_series = g_total / g_calor + plot_layout(heights = c(1, 2.2))
fig(g_series, "series.pdf", 6.5, 6)

L = log1p(Ytr)
R = cor(L)
dimnames(R) = list(colnames(Y), colnames(Y))
g_cor = ggplot(as.data.frame(as.table(R)), aes(factor(Var1, ordem), factor(Var2, ordem), fill = Freq)) + theme_bw(base_size = 9) +
  geom_tile() + scale_fill_gradient2(limits = c(-1, 1), low = "#2166ac", high = "#b2182b", name = "r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fig(g_cor, "cor.pdf", 4.3, 4)

#imagens: ACF e correlacao cruzada
sel = colnames(Y)[match(c(3509502, 3552403, 3520509, 3532009), cod)]
acf_df = bind_rows(lapply(sel, function(m) { a = acf(L[, m], lag.max = 104, plot = FALSE)
  data.frame(mun = m, lag = a$lag[-1], r = a$acf[-1]) }))
g_acf_s = ggplot(acf_df, aes(lag, r)) + theme_bw(base_size = 9) + geom_col(width = 0.4) +
  geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(Tn), linetype = 2, colour = "blue") +
  geom_vline(xintercept = 52, linetype = 3, colour = "red") + facet_wrap(~factor(mun, sel), ncol = 1) +
  labs(x = "Defasagem (semanas)", y = "FAC")
camp = which(cod == 3509502)
dL = diff(L)
ccf_df = bind_rows(lapply(setdiff(1:N, camp), function(i) {
  cc = ccf(dL[, camp], dL[, i], lag.max = 8, plot = FALSE)
  k = which.max(cc$acf)
  data.frame(mun = colnames(Y)[i], lag = -cc$lag[k], r = cc$acf[k]) }))
g_ccf = ggplot(ccf_df, aes(lag, reorder(mun, lag), fill = r)) + theme_bw(base_size = 9) + geom_col() + scale_fill_viridis_c(name = "r máx.") +
  labs(x = "Defasagem de máxima correlação com Campinas (semanas)", y = NULL)
g_acf_ccf = g_acf_s + g_ccf + plot_layout(widths = c(1, 1.2))
fig(g_acf_ccf, "acf_ccf.pdf", 6.5, 5)

#dados: estacionariedade (ADF e KPSS)
tab("Estacionariedade (p-valores)", data.frame(Municipio = colnames(Y),
  ADF = sapply(1:N, function(i) suppressWarnings(adf.test(L[, i])$p.value)),
  KPSS = sapply(1:N, function(i) suppressWarnings(kpss.test(L[, i])$p.value)))[order(-pop), ])

#imagens: media x variancia
bloco = rep(seq_len(ceiling(Tn / 13)), each = 13)[1:Tn]
mv = bind_rows(lapply(1:N, function(i) data.frame(m = tapply(Ytr[, i], bloco, mean), v = tapply(Ytr[, i], bloco, var)))) %>%
  filter(m > 0, v > 0)
cf = coef(lm(log(v) ~ log(m), data = mv))
incl = cf[2]
cat(sprintf("Inclinação média-variância = %.2f\n", incl))
ref = data.frame(m = range(mv$m))
g_media_variancia = ggplot(mv, aes(m, v)) + theme_bw(base_size = 9) + geom_point(size = 0.6, alpha = 0.4) +
  geom_line(data = rbind(transform(ref, v = m, curva = "Poisson: inclinação 1"),
                         transform(ref, v = m^2, curva = "Binomial negativa: inclinação 2"),
                         transform(ref, v = exp(cf[1]) * m^incl, curva = sprintf("Ajustada: inclinação %.2f", incl))),
            aes(m, v, colour = curva, linetype = curva), linewidth = 0.7) +
  scale_colour_manual(values = c("#b2182b", "black", "#2166ac"), name = NULL) +
  scale_linetype_manual(values = c(2, 1, 2), name = NULL) + scale_x_log10() + scale_y_log10() +
  labs(x = "Média no bloco de 13 semanas", y = "Variância no bloco") +
  theme(legend.position = "bottom", legend.direction = "vertical")
fig(g_media_variancia, "media_variancia.pdf", 3.8, 4.2)

#estimacao: ordem, rede e coeficientes
aj = ajusta_grade(r)
grade = aj$grade
melhor = aj$melhor
fit_best = aj$best
tab("Modelos escolhidos", melhor)
#imagens: criterios de selecao
g_criterios = grade %>% pivot_longer(c(AIC, BIC, QIC), names_to = "criterio", values_to = "valor") %>%
  ggplot(aes(p, valor, colour = rede)) + theme_bw(base_size = 9) + geom_line() + geom_point() + facet_grid(criterio ~ tipo, scales = "free_y") +
  labs(x = "Ordem p", y = NULL, colour = NULL) + theme(legend.position = "bottom")
fig(g_criterios, "criterios.pdf", 6.5, 6)
#dados: comparacao de redes e de ordens (diferenca para o minimo)
tab("Redes (na ordem escolhida)", grade %>% semi_join(melhor %>% select(tipo, p), by = c("tipo", "p")) %>%
  group_by(tipo) %>% mutate(dAIC = AIC - min(AIC), dQIC = QIC - min(QIC)) %>% ungroup() %>% select(tipo, p, rede, dAIC, dQIC), 1)
tab("Ordem p", grade %>% group_by(tipo, p) %>% summarise(AIC = min(AIC), BIC = min(BIC), QIC = min(QIC), .groups = "drop") %>%
  group_by(tipo) %>% mutate(dAIC = AIC - min(AIC), dBIC = BIC - min(BIC), dQIC = QIC - min(QIC)) %>% ungroup() %>%
  select(tipo, p, dAIC, dBIC, dQIC), 1)
tab_coef = tabela_coef(fit_best)
tab("Coeficientes", tab_coef, 4)
for (m in names(fit_best)) {
  b = fit_best[[m]]$fit$coefs[, 1]
  p = fit_best[[m]]$p
  cat(sprintf("%s: soma dos coeficientes = %.3f\n", m, if (m == "Linear") sum(b[2:(2*p+1)]) else sum(abs(b[2:(2*p+1)])))) }

#extensao sazonal
W_ll = fit_best[["Log-linear"]]$W
p_ll = fit_best[["Log-linear"]]$p
saz = lapply(0:2, function(K) pnar_sazonal(r, W_ll, p_ll, K))
tab("Sazonal", data.frame(Modelo = c("Pacote PNAR", "Propria K=0", "K=1", "K=2"),
  soma = c(sum(abs(fit_best[["Log-linear"]]$fit$coefs[2:(2*p_ll+1), 1])), sapply(saz, `[[`, "soma")),
  QIC = c(fit_best[["Log-linear"]]$fit$ic["QIC"], sapply(saz, `[[`, "QIC"))))
K_best = which.min(sapply(saz[2:3], `[[`, "QIC"))
tab("Coeficientes sazonais", cbind(par = rownames(saz[[K_best + 1]]$coefs), saz[[K_best + 1]]$coefs), 4)

#diagnostico dos residuos
diag_res = lapply(setNames(names(fit_best), names(fit_best)), function(m) {
  idx = (fit_best[[m]]$p + 1):Tn
  lam = lambda_pnar(fit_best[[m]], m, Y, idx)
  (Y[idx, ] - lam) / sqrt(lam) })
tab("Diagnostico", bind_rows(lapply(names(diag_res), function(m) { E = diag_res[[m]]
  data.frame(Modelo = m, Dispersao = sum(E^2) / (length(E) - nrow(fit_best[[m]]$fit$coefs)),
    LB_e = sum(apply(E, 2, function(x) Box.test(x, 10, "Ljung-Box")$p.value < 0.05)),
    LB_e2 = sum(apply(E, 2, function(x) Box.test(x^2, 10, "Ljung-Box")$p.value < 0.05))) })))
#imagens: residuos
Em = diag_res[["Log-linear"]]
acf_res = sapply(1:26, function(l) mean(sapply(1:N, function(i) cor(Em[-(1:l), i], Em[1:(nrow(Em) - l), i]))))
Rres = cor(Em)
dimnames(Rres) = list(colnames(Y), colnames(Y))
g_diagnostico = ggplot(data.frame(lag = 1:26, r = acf_res), aes(lag, r)) + theme_bw(base_size = 9) + geom_col(width = 0.3) +
    geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(nrow(Em)), linetype = 2, colour = "blue") +
    labs(x = "Defasagem (semanas)", y = "ACF média dos resíduos") +
  ggplot(as.data.frame(as.table(Rres)), aes(factor(Var1, ordem), factor(Var2, ordem), fill = Freq)) + theme_bw(base_size = 9) + geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#2166ac", high = "#b2182b", name = "r") + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6)) +
  plot_layout(widths = c(1, 1.3))
fig(g_diagnostico, "diagnostico.pdf", 6.5, 3.3)
fl = fit_best[["Linear"]]
cat("\n## Teste de linearidade\n")
print(score_test_nonlinpq_h0(fl$fit$coefs[, 1], Ytr, fl$W, fl$p, d = 1, Z = fl$Z))

#previsao 2024
pv = preve(r, fit_best, list("PNAR log-linear sazonal" = prev_sazonal(saz[[K_best + 1]], r, which(!r$treino))))
tab("Previsao 2024 (RMC)", pv$tab, 2)
#imagens: previsao
idx_te = pv$idx
mods = c(paste("PNAR", names(fit_best)), "PNAR log-linear sazonal")
dp = bind_rows(lapply(mods, function(m) data.frame(data = rep(datas[idx_te], length(sel)),
  mun = rep(sel, each = length(idx_te)), mu = c(pv$prev[[m]][, match(sel, colnames(Y))]), modelo = m)))
do = data.frame(data = rep(datas[idx_te], length(sel)), mun = rep(sel, each = length(idx_te)), y = c(Y[idx_te, sel]))
g_previsao = ggplot() + theme_bw(base_size = 9) + geom_point(data = do, aes(data, y), size = 0.6) +
  geom_line(data = dp, aes(data, mu, colour = modelo), linewidth = 0.4) +
  facet_wrap(~factor(mun, sel), scales = "free_y") + scale_x_date(date_labels = "%b") +
  labs(x = NULL, y = "Casos semanais (2024)", colour = NULL) + theme(legend.position = "bottom")
fig(g_previsao, "previsao.pdf", 6.5, 4.5)

#sensibilidade: regiao ampliada (87 municipios)
g = prepara(cod_grande)
aj_g = ajusta_grade(g)
cat(sprintf("\nRegião maior: N = %d\n", g$N))
tab("Redes (regiao ampliada)", tabela_redes(g), 2)
tab("Modelos escolhidos (ampliada)", aj_g$melhor)
tab("Coeficientes (ampliada)", tabela_coef(aj_g$best), 4)
saz_g = pnar_sazonal(g, aj_g$best[["Log-linear"]]$W, aj_g$best[["Log-linear"]]$p, K_best)
pv_g = preve(g, aj_g$best, list("PNAR log-linear sazonal" = prev_sazonal(saz_g, g, which(!g$treino))))
tab("Previsao 2024 (ampliada)", pv_g$tab, 2)
#imagens: mapa da regiao ampliada
g_mapa_grande = ggplot(cbind(g$malha_plot, inc = 1e5 * colMeans(g$Ytr) / g$pop)) + theme_bw(base_size = 9) +
  geom_sf(aes(fill = inc), colour = "grey40", linewidth = 0.1) +
  geom_sf(data = st_union(malha_plot[malha_plot$cod %in% cod_rmc, ]), fill = NA, colour = "red", linewidth = 0.5) +
  scale_fill_viridis_c(trans = "log10", name = "Casos/100 mil\npor semana") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
fig(g_mapa_grande, "mapa_grande.pdf", 5.5, 4.2)
