function Hminus(B, B0, n)
    1 / (1 + (B / B0)^n)
end

function Hs(B, B0, n, fc)
    Hminus(B, B0, n) + fc * (1 - Hminus(B, B0, n))
end

function model(du, u, p, t)
    k_200, k_ZEB, k_SNAI, k_SLUG, k_CDH, d_200, dR_ZEB, dR_SNAI, dR_SLUG, dR_CDH, l_ZEB, l_SNAI, l_SLUG, l_CDH, dP_ZEB, dP_SNAI, dP_SLUG, dP_CDH, K_200_ZEB, K_200_SNAI, K_200_SLUG, n_200_ZEB, n_200_SNAI, n_200_SLUG, fc_200_ZEB, fc_200_SNAI, fc_200_SLUG, K_ZEB_ZEB, K_ZEB_SNAI, K_ZEB_SLUG, n_ZEB_ZEB, n_ZEB_SNAI, n_ZEB_SLUG, fc_ZEB_ZEB, fc_ZEB_SNAI, fc_ZEB_SLUG, K_SNAI_SNAI, K_SNAI_SLUG, n_SNAI_SNAI, n_SNAI_SLUG, fc_SNAI_SNAI, fc_SNAI_SLUG, K_SLUG_SLUG, K_SLUG_SNAI, n_SLUG_SLUG, n_SLUG_SNAI, fc_SLUG_SLUG, fc_SLUG_SNAI, K_CDH_SLUG, K_CDH_SNAI, K_CDH_ZEB, n_CDH_SLUG, n_CDH_SNAI, n_CDH_ZEB, fc_CDH_SLUG, fc_CDH_SNAI, fc_CDH_ZEB, t_ZEB_1, t_ZEB_2, t_SLUG, a_ZEB_1, a_ZEB_2, a_SLUG, b_ZEB_1, b_ZEB_2, b_SLUG, kOff_ZEB, kOff_SLUG, KA_ZEB, KA_SLUG, σR_200, σR_ZEB, σR_SNAI, σR_SLUG, σR_CDH, σX_ZEB, σX_SNAI, σX_SLUG, σX_CDH, σC1_ZEB, σC2_ZEB, σC1_SLUG, τ = p
    R_200, R_ZEB, R_SNAI, R_SLUG, R_CDH, X_ZEB, X_SNAI, X_SLUG, X_CDH, C1_ZEB, C2_ZEB, C1_SLUG, ξR_200, ξR_ZEB, ξR_SNAI, ξR_SLUG, ξR_CDH, ξX_ZEB, ξX_SNAI, ξX_SLUG, ξX_CDH, ξC1_ZEB, ξC2_ZEB, ξC1_SLUG = u
    R_200 = max(R_200, 0.)
    R_ZEB = max(R_ZEB, 0.)
    R_SNAI = max(R_SNAI, 0.)
    R_SLUG = max(R_SLUG, 0.)
    R_CDH = max(R_CDH, 0.)
    X_ZEB = max(X_ZEB, 0.)
    X_SNAI = max(X_SNAI, 0.)
    X_SLUG = max(X_SLUG, 0.)
    X_CDH = max(X_CDH, 0.)
    C1_ZEB = max(C1_ZEB, 0.)
    C2_ZEB = max(C2_ZEB, 0.)
    C1_SLUG = max(C1_SLUG, 0.)

    kOn_ZEB = kOff_ZEB * KA_ZEB
    kOn_SLUG = kOff_SLUG * KA_SLUG

    JI = k_200 * Hs(X_ZEB, K_200_ZEB, n_200_ZEB, fc_200_ZEB) * Hs(X_SNAI, K_200_SNAI, n_200_SNAI, fc_200_SNAI) * Hs(X_SLUG, K_200_SLUG, n_200_SLUG, fc_200_SLUG) - d_200 * R_200
    JRZ = k_ZEB * Hs(X_ZEB, K_ZEB_ZEB, n_ZEB_ZEB, fc_ZEB_ZEB) * Hs(X_SNAI, K_ZEB_SNAI, n_ZEB_SNAI, fc_ZEB_SNAI) * Hs(X_SLUG, K_ZEB_SLUG, n_ZEB_SLUG, fc_ZEB_SLUG) - dR_ZEB * R_ZEB
    JRSn = k_SNAI * Hs(X_SNAI, K_SNAI_SNAI, n_SNAI_SNAI, fc_SNAI_SNAI) * Hs(X_SLUG, K_SNAI_SLUG, n_SNAI_SLUG, fc_SNAI_SLUG) - dR_SNAI * R_SNAI
    JRSl = k_SLUG * Hs(X_SLUG, K_SLUG_SLUG, n_SLUG_SLUG, fc_SLUG_SLUG) * Hs(X_SNAI, K_SLUG_SNAI, n_SLUG_SNAI, fc_SLUG_SNAI) - dR_SLUG * R_SLUG
    JRC = k_CDH * Hs(X_SLUG, K_CDH_SLUG, n_CDH_SLUG, fc_CDH_SLUG) * Hs(X_SNAI, K_CDH_SNAI, n_CDH_SNAI, fc_CDH_SNAI) * Hs(X_ZEB, K_CDH_ZEB, n_CDH_ZEB, fc_CDH_ZEB) - dR_CDH * R_CDH
    JPZ = l_ZEB * (R_ZEB + t_ZEB_1 * C1_ZEB + t_ZEB_2 * C2_ZEB) - dP_ZEB * X_ZEB
    JPSn = l_SNAI * R_SNAI - dP_SNAI * X_SNAI
    JPSl = l_SLUG * (R_SLUG + t_SLUG * C1_SLUG) - dP_SLUG * X_SLUG
    JPC = l_CDH * R_CDH - dP_CDH * X_CDH
    JC1Z = 2 * kOn_ZEB * R_ZEB * R_200 - kOff_ZEB * C1_ZEB
    JC1ZDR = dR_ZEB * a_ZEB_1 * C1_ZEB
    JC1ZDI = d_200 * b_ZEB_1 * C1_ZEB
    JC2Z = kOn_ZEB * C1_ZEB * R_200 - 2 * kOff_ZEB * C2_ZEB
    JC2ZDR = dR_ZEB * a_ZEB_2 * C2_ZEB
    JC2ZDI = 2 * d_200 * b_ZEB_2 * C2_ZEB
    JCSl = kOn_SLUG * R_SLUG * R_200 - kOff_SLUG * C1_SLUG
    JCSlDR = dR_SLUG * a_SLUG * C1_SLUG
    JCSlDI = d_200 * b_SLUG * C1_SLUG

    du[1] = (1)*JI + (-1)*JC1Z + (1)*JC1ZDR + (-1)*JC2Z + (2)*JC2ZDR + (-1)*JCSl + (1)*JCSlDR + ξR_200
    du[2] = (1)*JRZ + (-1)*JC1Z + (1)*JC1ZDI + ξR_ZEB
    du[3] = (1)*JRSn + ξR_SNAI
    du[4] = (1)*JRSl + (-1)*JCSl + (1)*JCSlDI + ξR_SLUG
    du[5] = (1)*JRC + ξR_CDH
    du[6] = (1)*JPZ + ξX_ZEB
    du[7] = (1)*JPSn + ξX_SNAI
    du[8] = (1)*JPSl + ξX_SLUG
    du[9] = (1)*JPC + ξX_CDH
    du[10] = (1)*JC1Z + (-1)*JC1ZDR + (-1)*JC1ZDI + (-1)*JC2Z + (1)*JC2ZDI + ξC1_ZEB
    du[11] = (1)*JC2Z + (-1)*JC2ZDR + (-1)*JC2ZDI + ξC2_ZEB
    du[12] = (1)*JCSl + (-1)*JCSlDR + (-1)*JCSlDI + ξC1_SLUG
    du[13] = -u[13] / τ
    du[14] = -u[14] / τ
    du[15] = -u[15] / τ
    du[16] = -u[16] / τ
    du[17] = -u[17] / τ
    du[18] = -u[18] / τ
    du[19] = -u[19] / τ
    du[20] = -u[20] / τ
    du[21] = -u[21] / τ
    du[22] = -u[22] / τ
    du[23] = -u[23] / τ
    du[24] = -u[24] / τ
end

function stoch(du, u, p, t)
    R_200, R_ZEB, R_SNAI, R_SLUG, R_CDH, X_ZEB, X_SNAI, X_SLUG, X_CDH, C1_ZEB, C2_ZEB, C1_SLUG, ξR_200, ξR_ZEB, ξR_SNAI, ξR_SLUG, ξR_CDH, ξX_ZEB, ξX_SNAI, ξX_SLUG, ξX_CDH, ξC1_ZEB, ξC2_ZEB, ξC1_SLUG = u
    k_200, k_ZEB, k_SNAI, k_SLUG, k_CDH, d_200, dR_ZEB, dR_SNAI, dR_SLUG, dR_CDH, l_ZEB, l_SNAI, l_SLUG, l_CDH, dP_ZEB, dP_SNAI, dP_SLUG, dP_CDH, K_200_ZEB, K_200_SNAI, K_200_SLUG, n_200_ZEB, n_200_SNAI, n_200_SLUG, fc_200_ZEB, fc_200_SNAI, fc_200_SLUG, K_ZEB_ZEB, K_ZEB_SNAI, K_ZEB_SLUG, n_ZEB_ZEB, n_ZEB_SNAI, n_ZEB_SLUG, fc_ZEB_ZEB, fc_ZEB_SNAI, fc_ZEB_SLUG, K_SNAI_SNAI, K_SNAI_SLUG, n_SNAI_SNAI, n_SNAI_SLUG, fc_SNAI_SNAI, fc_SNAI_SLUG, K_SLUG_SLUG, K_SLUG_SNAI, n_SLUG_SLUG, n_SLUG_SNAI, fc_SLUG_SLUG, fc_SLUG_SNAI, K_CDH_SLUG, K_CDH_SNAI, K_CDH_ZEB, n_CDH_SLUG, n_CDH_SNAI, n_CDH_ZEB, fc_CDH_SLUG, fc_CDH_SNAI, fc_CDH_ZEB, t_ZEB_1, t_ZEB_2, t_SLUG, a_ZEB_1, a_ZEB_2, a_SLUG, b_ZEB_1, b_ZEB_2, b_SLUG, kOff_ZEB, kOff_SLUG, KA_ZEB, KA_SLUG, σR_200, σR_ZEB, σR_SNAI, σR_SLUG, σR_CDH, σX_ZEB, σX_SNAI, σX_SLUG, σX_CDH, σC1_ZEB, σC2_ZEB, σC1_SLUG, τ = p

    du[1] = 0
    du[2] = 0
    du[3] = 0
    du[4] = 0
    du[5] = 0
    du[6] = 0
    du[7] = 0
    du[8] = 0
    du[9] = 0
    du[10] = 0
    du[11] = 0
    du[12] = 0
    du[13] = σR_200 * sqrt(2.0 / τ) * R_200
    du[14] = σR_ZEB * sqrt(2.0 / τ) * R_ZEB
    du[15] = σR_SNAI * sqrt(2.0 / τ) * R_SNAI
    du[16] = σR_SLUG * sqrt(2.0 / τ) * R_SLUG
    du[17] = σR_CDH * sqrt(2.0 / τ) * R_CDH
    du[18] = σX_ZEB * sqrt(2.0 / τ) * X_ZEB
    du[19] = σX_SNAI * sqrt(2.0 / τ) * X_SNAI
    du[20] = σX_SLUG * sqrt(2.0 / τ) * X_SLUG
    du[21] = σX_CDH * sqrt(2.0 / τ) * X_CDH
    du[22] = σC1_ZEB * sqrt(2.0 / τ) * C1_ZEB
    du[23] = σC2_ZEB * sqrt(2.0 / τ) * C2_ZEB
    du[24] = σC1_SLUG * sqrt(2.0 / τ) * C1_SLUG
end

species_ids = Dict(
    "R_200" => 1,
    "R_ZEB" => 2,
    "R_SNAI" => 3,
    "R_SLUG" => 4,
    "R_CDH" => 5,
    "X_ZEB" => 6,
    "X_SNAI" => 7,
    "X_SLUG" => 8,
    "X_CDH" => 9,
    "C1_ZEB" => 10,
    "C2_ZEB" => 11,
    "C1_SLUG" => 12,
    "ξR_200" => 13,
    "ξR_ZEB" => 14,
    "ξR_SNAI" => 15,
    "ξR_SLUG" => 16,
    "ξR_CDH" => 17,
    "ξX_ZEB" => 18,
    "ξX_SNAI" => 19,
    "ξX_SLUG" => 20,
    "ξX_CDH" => 21,
    "ξC1_ZEB" => 22,
    "ξC2_ZEB" => 23,
    "ξC1_SLUG" => 24,
)

param_ids = Dict(
    "k_200" => 1,
    "k_ZEB" => 2,
    "k_SNAI" => 3,
    "k_SLUG" => 4,
    "k_CDH" => 5,
    "d_200" => 6,
    "dR_ZEB" => 7,
    "dR_SNAI" => 8,
    "dR_SLUG" => 9,
    "dR_CDH" => 10,
    "l_ZEB" => 11,
    "l_SNAI" => 12,
    "l_SLUG" => 13,
    "l_CDH" => 14,
    "dP_ZEB" => 15,
    "dP_SNAI" => 16,
    "dP_SLUG" => 17,
    "dP_CDH" => 18,
    "K_200_ZEB" => 19,
    "K_200_SNAI" => 20,
    "K_200_SLUG" => 21,
    "n_200_ZEB" => 22,
    "n_200_SNAI" => 23,
    "n_200_SLUG" => 24,
    "fc_200_ZEB" => 25,
    "fc_200_SNAI" => 26,
    "fc_200_SLUG" => 27,
    "K_ZEB_ZEB" => 28,
    "K_ZEB_SNAI" => 29,
    "K_ZEB_SLUG" => 30,
    "n_ZEB_ZEB" => 31,
    "n_ZEB_SNAI" => 32,
    "n_ZEB_SLUG" => 33,
    "fc_ZEB_ZEB" => 34,
    "fc_ZEB_SNAI" => 35,
    "fc_ZEB_SLUG" => 36,
    "K_SNAI_SNAI" => 37,
    "K_SNAI_SLUG" => 38,
    "n_SNAI_SNAI" => 39,
    "n_SNAI_SLUG" => 40,
    "fc_SNAI_SNAI" => 41,
    "fc_SNAI_SLUG" => 42,
    "K_SLUG_SLUG" => 43,
    "K_SLUG_SNAI" => 44,
    "n_SLUG_SLUG" => 45,
    "n_SLUG_SNAI" => 46,
    "fc_SLUG_SLUG" => 47,
    "fc_SLUG_SNAI" => 48,
    "K_CDH_SLUG" => 49,
    "K_CDH_SNAI" => 50,
    "K_CDH_ZEB" => 51,
    "n_CDH_SLUG" => 52,
    "n_CDH_SNAI" => 53,
    "n_CDH_ZEB" => 54,
    "fc_CDH_SLUG" => 55,
    "fc_CDH_SNAI" => 56,
    "fc_CDH_ZEB" => 57,
    "t_ZEB_1" => 58,
    "t_ZEB_2" => 59,
    "t_SLUG" => 60,
    "a_ZEB_1" => 61,
    "a_ZEB_2" => 62,
    "a_SLUG" => 63,
    "b_ZEB_1" => 64,
    "b_ZEB_2" => 65,
    "b_SLUG" => 66,
    "kOff_ZEB" => 67,
    "kOff_SLUG" => 68,
    "KA_ZEB" => 69,
    "KA_SLUG" => 70,
    "σR_200" => 71,
    "σR_ZEB" => 72,
    "σR_SNAI" => 73,
    "σR_SLUG" => 74,
    "σR_CDH" => 75,
    "σX_ZEB" => 76,
    "σX_SNAI" => 77,
    "σX_SLUG" => 78,
    "σX_CDH" => 79,
    "σC1_ZEB" => 80,
    "σC2_ZEB" => 81,
    "σC1_SLUG" => 82,
    "τ" => 83,
)

p_default = [
    10, # 1 k_200
    0.01, # 2 k_ZEB
    0.1, # 3 k_SNAI
    1, # 4 k_SLUG
    10, # 5 k_CDH
    1, # 6 d_200
    1, # 7 dR_ZEB
    1, # 8 dR_SNAI
    1, # 9 dR_SLUG
    1, # 10 dR_CDH
    1, # 11 l_ZEB
    1.8, # 12 l_SNAI
    5, # 13 l_SLUG
    1, # 14 l_CDH
    1, # 15 dP_ZEB
    1.25, # 16 dP_SNAI
    1.1, # 17 dP_SLUG
    1, # 18 dP_CDH
    2.2, # 19 K_200_ZEB
    1.9, # 20 K_200_SNAI
    2.2, # 21 K_200_SLUG
    3, # 22 n_200_ZEB
    2, # 23 n_200_SNAI
    1, # 24 n_200_SLUG
    0.1, # 25 fc_200_ZEB
    0.1, # 26 fc_200_SNAI
    0.4, # 27 fc_200_SLUG
    0.25, # 28 K_ZEB_ZEB
    1.8, # 29 K_ZEB_SNAI
    2, # 30 K_ZEB_SLUG
    2, # 31 n_ZEB_ZEB
    2, # 32 n_ZEB_SNAI
    2, # 33 n_ZEB_SLUG
    7.5, # 34 fc_ZEB_ZEB
    10, # 35 fc_ZEB_SNAI
    4, # 36 fc_ZEB_SLUG
    1.8, # 37 K_SNAI_SNAI
    2.25, # 38 K_SNAI_SLUG
    5, # 39 n_SNAI_SNAI
    3, # 40 n_SNAI_SLUG
    0.4, # 41 fc_SNAI_SNAI
    0.5, # 42 fc_SNAI_SLUG
    2.5, # 43 K_SLUG_SLUG
    1.8, # 44 K_SLUG_SNAI
    4, # 45 n_SLUG_SLUG
    1, # 46 n_SLUG_SNAI
    4, # 47 fc_SLUG_SLUG
    0.5, # 48 fc_SLUG_SNAI
    0.5, # 49 K_CDH_SLUG
    0.5, # 50 K_CDH_SNAI
    1, # 51 K_CDH_ZEB
    2, # 52 n_CDH_SLUG
    2, # 53 n_CDH_SNAI
    2, # 54 n_CDH_ZEB
    0.3, # 55 fc_CDH_SLUG
    0.3, # 56 fc_CDH_SNAI
    0.75, # 57 fc_CDH_ZEB
    0.1, # 58 t_ZEB_1
    0.1, # 59 t_ZEB_2
    0.1, # 60 t_SLUG
    2, # 61 a_ZEB_1
    4, # 62 a_ZEB_2
    2, # 63 a_SLUG
    0.5, # 64 b_ZEB_1
    0.25, # 65 b_ZEB_2
    0.5, # 66 b_SLUG
    100, # 67 kOff_ZEB
    100, # 68 kOff_SLUG
    10, # 69 KA_ZEB
    10, # 70 KA_SLUG
    0, # 71 σR_200
    0, # 72 σR_ZEB
    0, # 73 σR_SNAI
    0, # 74 σR_SLUG
    0, # 75 σR_CDH
    0, # 76 σX_ZEB
    0, # 77 σX_SNAI
    0, # 78 σX_SLUG
    0, # 79 σX_CDH
    0, # 80 σC1_ZEB
    0, # 81 σC2_ZEB
    0, # 82 σC1_SLUG
    0.1, # 83 τ
]
