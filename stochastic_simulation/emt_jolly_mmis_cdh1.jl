function Hminus(B, B0, n)
    1 / (1 + (B / B0)^n)
end

function Hs(B, B0, n, fc)
    Hminus(B, B0, n) + fc * (1 - Hminus(B, B0, n))
end

function model(du, u, p, t)
    k_101, k_200, k_ZEB, k_SNAI, k_SLUG, k_CDH, d_101, d_200, dR_ZEB, dR_SNAI, dR_SLUG, dR_CDH, l_ZEB, l_SNAI, l_SLUG, l_CDH, dP_ZEB, dP_SNAI, dP_SLUG, dP_CDH, K_101_SNAI, K_101_SLUG, n_101_SNAI, n_101_SLUG, fc_101_SNAI, fc_101_SLUG, K_200_ZEB, K_200_SNAI, K_200_SLUG, n_200_ZEB, n_200_SNAI, n_200_SLUG, fc_200_ZEB, fc_200_SNAI, fc_200_SLUG, K_ZEB_ZEB, K_ZEB_SNAI, K_ZEB_SLUG, n_ZEB_ZEB, n_ZEB_SNAI, n_ZEB_SLUG, fc_ZEB_ZEB, fc_ZEB_SNAI, fc_ZEB_SLUG, K_SNAI_SNAI, K_SNAI_SLUG, n_SNAI_SNAI, n_SNAI_SLUG, fc_SNAI_SNAI, fc_SNAI_SLUG, K_SLUG_SLUG, K_SLUG_SNAI, n_SLUG_SLUG, n_SLUG_SNAI, fc_SLUG_SLUG, fc_SLUG_SNAI, K_CDH_SLUG, K_CDH_SNAI, K_CDH_ZEB, n_CDH_SLUG, n_CDH_SNAI, n_CDH_ZEB, fc_CDH_SLUG, fc_CDH_SNAI, fc_CDH_ZEB, t_ZEB_2, t_ZEB_22, t_ZEB_1, t_ZEB_12, t_ZEB_122, t_SLUG, a_ZEB_2, a_ZEB_22, a_ZEB_1, a_ZEB_12, a_ZEB_122, a_SLUG, b_200_ZEB_2, b_200_ZEB_22, b_200_ZEB_12, b_200_ZEB_122, b_101_ZEB_1, b_101_ZEB_12, b_101_ZEB_122, b_SLUG, kOff, KA_ZEB_101, KA_ZEB_200, KA_SLUG, σR_101, σR_200, σR_ZEB, σR_SNAI, σR_SLUG, σR_CDH, σX_ZEB, σX_SNAI, σX_SLUG, σX_CDH, σC_ZEB_2, σC_ZEB_22, σC_ZEB_1, σC_ZEB_12, σC_ZEB_122, σC1_SLUG, τ = p
    R_101, R_200, R_ZEB, R_SNAI, R_SLUG, R_CDH, X_ZEB, X_SNAI, X_SLUG, X_CDH, C_ZEB_2, C_ZEB_22, C_ZEB_1, C_ZEB_12, C_ZEB_122, C1_SLUG, ξR_101, ξR_200, ξR_ZEB, ξR_SNAI, ξR_SLUG, ξR_CDH, ξX_ZEB, ξX_SNAI, ξX_SLUG, ξX_CDH, ξC_ZEB_2, ξC_ZEB_22, ξC_ZEB_1, ξC_ZEB_12, ξC_ZEB_122, ξC1_SLUG = u
    R_101 = max(R_101, 0.)
    R_200 = max(R_200, 0.)
    R_ZEB = max(R_ZEB, 0.)
    R_SNAI = max(R_SNAI, 0.)
    R_SLUG = max(R_SLUG, 0.)
    R_CDH = max(R_CDH, 0.)
    X_ZEB = max(X_ZEB, 0.)
    X_SNAI = max(X_SNAI, 0.)
    X_SLUG = max(X_SLUG, 0.)
    X_CDH = max(X_CDH, 0.)
    C_ZEB_2 = max(C_ZEB_2, 0.)
    C_ZEB_22 = max(C_ZEB_22, 0.)
    C_ZEB_1 = max(C_ZEB_1, 0.)
    C_ZEB_12 = max(C_ZEB_12, 0.)
    C_ZEB_122 = max(C_ZEB_122, 0.)
    C1_SLUG = max(C1_SLUG, 0.)

    kOn_ZEB_101 = kOff * KA_ZEB_101
    kOn_ZEB_200 = kOff * KA_ZEB_200
    kOn_SLUG = kOff * KA_SLUG

    JI101 = k_101 * Hs(X_SNAI, K_101_SNAI, n_101_SNAI, fc_101_SNAI) * Hs(X_SLUG, K_101_SLUG, n_101_SLUG, fc_101_SLUG) - d_101 * R_101
    JI200 = k_200 * Hs(X_ZEB, K_200_ZEB, n_200_ZEB, fc_200_ZEB) * Hs(X_SNAI, K_200_SNAI, n_200_SNAI, fc_200_SNAI) * Hs(X_SLUG, K_200_SLUG, n_200_SLUG, fc_200_SLUG) - d_200 * R_200
    JRZ = k_ZEB * Hs(X_ZEB, K_ZEB_ZEB, n_ZEB_ZEB, fc_ZEB_ZEB) * Hs(X_SNAI, K_ZEB_SNAI, n_ZEB_SNAI, fc_ZEB_SNAI) * Hs(X_SLUG, K_ZEB_SLUG, n_ZEB_SLUG, fc_ZEB_SLUG) - dR_ZEB * R_ZEB
    JRSn = k_SNAI * Hs(X_SNAI, K_SNAI_SNAI, n_SNAI_SNAI, fc_SNAI_SNAI) * Hs(X_SLUG, K_SNAI_SLUG, n_SNAI_SLUG, fc_SNAI_SLUG) - dR_SNAI * R_SNAI
    JRSl = k_SLUG * Hs(X_SLUG, K_SLUG_SLUG, n_SLUG_SLUG, fc_SLUG_SLUG) * Hs(X_SNAI, K_SLUG_SNAI, n_SLUG_SNAI, fc_SLUG_SNAI) - dR_SLUG * R_SLUG
    JRC = k_CDH * Hs(X_SLUG, K_CDH_SLUG, n_CDH_SLUG, fc_CDH_SLUG) * Hs(X_SNAI, K_CDH_SNAI, n_CDH_SNAI, fc_CDH_SNAI) * Hs(X_ZEB, K_CDH_ZEB, n_CDH_ZEB, fc_CDH_ZEB) - dR_CDH * R_CDH
    JPZ = l_ZEB * (R_ZEB + t_ZEB_2 * C_ZEB_2 + t_ZEB_22 * C_ZEB_22 + t_ZEB_1 * C_ZEB_1 + t_ZEB_12 * C_ZEB_12 + t_ZEB_122 * C_ZEB_122) - dP_ZEB * X_ZEB
    JPSn = l_SNAI * R_SNAI - dP_SNAI * X_SNAI
    JPSl = l_SLUG * (R_SLUG + t_SLUG * C1_SLUG) - dP_SLUG * X_SLUG
    JPC = l_CDH * R_CDH - dP_CDH * X_CDH
    JCZ2 = 2 * kOn_ZEB_200 * R_ZEB * R_200 - kOff * C_ZEB_2
    JCZ2DR = dR_ZEB * a_ZEB_2 * C_ZEB_2
    JCZ2DI = d_200 * b_200_ZEB_2 * C_ZEB_2
    JCZ22 = kOn_ZEB_200 * C_ZEB_2 * R_200 - 2 * kOff * C_ZEB_22
    JCZ22DR = dR_ZEB * a_ZEB_22 * C_ZEB_22
    JCZ22DI = 2 * d_200 * b_200_ZEB_22 * C_ZEB_22
    JCZ1 = kOn_ZEB_101 * R_ZEB * R_101 - kOff * C_ZEB_1
    JCZ1DR = dR_ZEB * a_ZEB_1 * C_ZEB_1
    JCZ1DI = d_101 * b_101_ZEB_1 * C_ZEB_1
    JCZ12 = 2 * kOn_ZEB_200 * C_ZEB_1 * R_200 - kOff * C_ZEB_12
    JCZ21 = kOn_ZEB_101 * C_ZEB_2 * R_101 - kOff * C_ZEB_12
    JCZ12DR = dR_ZEB * a_ZEB_12 * C_ZEB_12
    JCZ12DI1 = d_101 * b_101_ZEB_12 * C_ZEB_12
    JCZ12DI2 = d_200 * b_200_ZEB_12 * C_ZEB_12
    JCZ122 = kOn_ZEB_200 * C_ZEB_12 * R_200 - kOff * C_ZEB_122
    JCZ221 = kOn_ZEB_101 * C_ZEB_22 * R_101 - kOff * C_ZEB_122
    JCZ122DR = dR_ZEB * a_ZEB_122 * C_ZEB_122
    JCZ122DI1 = d_101 * b_101_ZEB_122 * C_ZEB_122
    JCZ122DI2 = 2 * d_200 * b_200_ZEB_122 * C_ZEB_122
    JCSl = kOn_SLUG * R_SLUG * R_200 - kOff * C1_SLUG
    JCSlDR = dR_SLUG * a_SLUG * C1_SLUG
    JCSlDI = d_200 * b_SLUG * C1_SLUG

    du[1] = (1)*JI101 + (-1)*JCZ1 + (1)*JCZ1DR + (-1)*JCZ21 + (1)*JCZ12DR + (-1)*JCZ221 + (1)*JCZ122DR + ξR_101
    du[2] = (1)*JI200 + (-1)*JCZ2 + (1)*JCZ2DR + (-1)*JCZ22 + (2)*JCZ22DR + (-1)*JCZ12 + (1)*JCZ12DR + (-1)*JCZ122 + (2)*JCZ122DR + (-1)*JCSl + (1)*JCSlDR + ξR_200
    du[3] = (1)*JRZ + (-1)*JCZ2 + (1)*JCZ2DI + (-1)*JCZ1 + (1)*JCZ1DI + ξR_ZEB
    du[4] = (1)*JRSn + ξR_SNAI
    du[5] = (1)*JRSl + (-1)*JCSl + (1)*JCSlDI + ξR_SLUG
    du[6] = (1)*JRC + ξR_CDH
    du[7] = (1)*JPZ + ξX_ZEB
    du[8] = (1)*JPSn + ξX_SNAI
    du[9] = (1)*JPSl + ξX_SLUG
    du[10] = (1)*JPC + ξX_CDH
    du[11] = (1)*JCZ2 + (-1)*JCZ2DR + (-1)*JCZ2DI + (-1)*JCZ22 + (1)*JCZ22DI + (-1)*JCZ21 + (1)*JCZ12DI1 + ξC_ZEB_2
    du[12] = (1)*JCZ22 + (-1)*JCZ22DR + (-1)*JCZ22DI + (-1)*JCZ221 + (1)*JCZ122DI1 + ξC_ZEB_22
    du[13] = (1)*JCZ1 + (-1)*JCZ1DR + (-1)*JCZ1DI + (-1)*JCZ12 + (1)*JCZ12DI2 + ξC_ZEB_1
    du[14] = (1)*JCZ12 + (1)*JCZ21 + (-1)*JCZ12DR + (-1)*JCZ12DI1 + (-1)*JCZ12DI2 + (-1)*JCZ122 + (1)*JCZ122DI2 + ξC_ZEB_12
    du[15] = (1)*JCZ122 + (1)*JCZ221 + (-1)*JCZ122DR + (-1)*JCZ122DI1 + (-1)*JCZ122DI2 + ξC_ZEB_122
    du[16] = (1)*JCSl + (-1)*JCSlDR + (-1)*JCSlDI + ξC1_SLUG
    du[17] = -u[17] / τ
    du[18] = -u[18] / τ
    du[19] = -u[19] / τ
    du[20] = -u[20] / τ
    du[21] = -u[21] / τ
    du[22] = -u[22] / τ
    du[23] = -u[23] / τ
    du[24] = -u[24] / τ
    du[25] = -u[25] / τ
    du[26] = -u[26] / τ
    du[27] = -u[27] / τ
    du[28] = -u[28] / τ
    du[29] = -u[29] / τ
    du[30] = -u[30] / τ
    du[31] = -u[31] / τ
    du[32] = -u[32] / τ
end

function stoch(du, u, p, t)
    R_101, R_200, R_ZEB, R_SNAI, R_SLUG, R_CDH, X_ZEB, X_SNAI, X_SLUG, X_CDH, C_ZEB_2, C_ZEB_22, C_ZEB_1, C_ZEB_12, C_ZEB_122, C1_SLUG, ξR_101, ξR_200, ξR_ZEB, ξR_SNAI, ξR_SLUG, ξR_CDH, ξX_ZEB, ξX_SNAI, ξX_SLUG, ξX_CDH, ξC_ZEB_2, ξC_ZEB_22, ξC_ZEB_1, ξC_ZEB_12, ξC_ZEB_122, ξC1_SLUG = u
    k_101, k_200, k_ZEB, k_SNAI, k_SLUG, k_CDH, d_101, d_200, dR_ZEB, dR_SNAI, dR_SLUG, dR_CDH, l_ZEB, l_SNAI, l_SLUG, l_CDH, dP_ZEB, dP_SNAI, dP_SLUG, dP_CDH, K_101_SNAI, K_101_SLUG, n_101_SNAI, n_101_SLUG, fc_101_SNAI, fc_101_SLUG, K_200_ZEB, K_200_SNAI, K_200_SLUG, n_200_ZEB, n_200_SNAI, n_200_SLUG, fc_200_ZEB, fc_200_SNAI, fc_200_SLUG, K_ZEB_ZEB, K_ZEB_SNAI, K_ZEB_SLUG, n_ZEB_ZEB, n_ZEB_SNAI, n_ZEB_SLUG, fc_ZEB_ZEB, fc_ZEB_SNAI, fc_ZEB_SLUG, K_SNAI_SNAI, K_SNAI_SLUG, n_SNAI_SNAI, n_SNAI_SLUG, fc_SNAI_SNAI, fc_SNAI_SLUG, K_SLUG_SLUG, K_SLUG_SNAI, n_SLUG_SLUG, n_SLUG_SNAI, fc_SLUG_SLUG, fc_SLUG_SNAI, K_CDH_SLUG, K_CDH_SNAI, K_CDH_ZEB, n_CDH_SLUG, n_CDH_SNAI, n_CDH_ZEB, fc_CDH_SLUG, fc_CDH_SNAI, fc_CDH_ZEB, t_ZEB_2, t_ZEB_22, t_ZEB_1, t_ZEB_12, t_ZEB_122, t_SLUG, a_ZEB_2, a_ZEB_22, a_ZEB_1, a_ZEB_12, a_ZEB_122, a_SLUG, b_200_ZEB_2, b_200_ZEB_22, b_200_ZEB_12, b_200_ZEB_122, b_101_ZEB_1, b_101_ZEB_12, b_101_ZEB_122, b_SLUG, kOff, KA_ZEB_101, KA_ZEB_200, KA_SLUG, σR_101, σR_200, σR_ZEB, σR_SNAI, σR_SLUG, σR_CDH, σX_ZEB, σX_SNAI, σX_SLUG, σX_CDH, σC_ZEB_2, σC_ZEB_22, σC_ZEB_1, σC_ZEB_12, σC_ZEB_122, σC1_SLUG, τ = p

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
    du[13] = 0
    du[14] = 0
    du[15] = 0
    du[16] = 0
    du[17] = σR_101 * sqrt(2.0 / τ) * R_101
    du[18] = σR_200 * sqrt(2.0 / τ) * R_200
    du[19] = σR_ZEB * sqrt(2.0 / τ) * R_ZEB
    du[20] = σR_SNAI * sqrt(2.0 / τ) * R_SNAI
    du[21] = σR_SLUG * sqrt(2.0 / τ) * R_SLUG
    du[22] = σR_CDH * sqrt(2.0 / τ) * R_CDH
    du[23] = σX_ZEB * sqrt(2.0 / τ) * X_ZEB
    du[24] = σX_SNAI * sqrt(2.0 / τ) * X_SNAI
    du[25] = σX_SLUG * sqrt(2.0 / τ) * X_SLUG
    du[26] = σX_CDH * sqrt(2.0 / τ) * X_CDH
    du[27] = σC_ZEB_2 * sqrt(2.0 / τ) * C_ZEB_2
    du[28] = σC_ZEB_22 * sqrt(2.0 / τ) * C_ZEB_22
    du[29] = σC_ZEB_1 * sqrt(2.0 / τ) * C_ZEB_1
    du[30] = σC_ZEB_12 * sqrt(2.0 / τ) * C_ZEB_12
    du[31] = σC_ZEB_122 * sqrt(2.0 / τ) * C_ZEB_122
    du[32] = σC1_SLUG * sqrt(2.0 / τ) * C1_SLUG
end

species_ids = Dict(
    "R_101" => 1,
    "R_200" => 2,
    "R_ZEB" => 3,
    "R_SNAI" => 4,
    "R_SLUG" => 5,
    "R_CDH" => 6,
    "X_ZEB" => 7,
    "X_SNAI" => 8,
    "X_SLUG" => 9,
    "X_CDH" => 10,
    "C_ZEB_2" => 11,
    "C_ZEB_22" => 12,
    "C_ZEB_1" => 13,
    "C_ZEB_12" => 14,
    "C_ZEB_122" => 15,
    "C1_SLUG" => 16,
    "ξR_101" => 17,
    "ξR_200" => 18,
    "ξR_ZEB" => 19,
    "ξR_SNAI" => 20,
    "ξR_SLUG" => 21,
    "ξR_CDH" => 22,
    "ξX_ZEB" => 23,
    "ξX_SNAI" => 24,
    "ξX_SLUG" => 25,
    "ξX_CDH" => 26,
    "ξC_ZEB_2" => 27,
    "ξC_ZEB_22" => 28,
    "ξC_ZEB_1" => 29,
    "ξC_ZEB_12" => 30,
    "ξC_ZEB_122" => 31,
    "ξC1_SLUG" => 32,
)

param_ids = Dict(
    "k_101" => 1,
    "k_200" => 2,
    "k_ZEB" => 3,
    "k_SNAI" => 4,
    "k_SLUG" => 5,
    "k_CDH" => 6,
    "d_101" => 7,
    "d_200" => 8,
    "dR_ZEB" => 9,
    "dR_SNAI" => 10,
    "dR_SLUG" => 11,
    "dR_CDH" => 12,
    "l_ZEB" => 13,
    "l_SNAI" => 14,
    "l_SLUG" => 15,
    "l_CDH" => 16,
    "dP_ZEB" => 17,
    "dP_SNAI" => 18,
    "dP_SLUG" => 19,
    "dP_CDH" => 20,
    "K_101_SNAI" => 21,
    "K_101_SLUG" => 22,
    "n_101_SNAI" => 23,
    "n_101_SLUG" => 24,
    "fc_101_SNAI" => 25,
    "fc_101_SLUG" => 26,
    "K_200_ZEB" => 27,
    "K_200_SNAI" => 28,
    "K_200_SLUG" => 29,
    "n_200_ZEB" => 30,
    "n_200_SNAI" => 31,
    "n_200_SLUG" => 32,
    "fc_200_ZEB" => 33,
    "fc_200_SNAI" => 34,
    "fc_200_SLUG" => 35,
    "K_ZEB_ZEB" => 36,
    "K_ZEB_SNAI" => 37,
    "K_ZEB_SLUG" => 38,
    "n_ZEB_ZEB" => 39,
    "n_ZEB_SNAI" => 40,
    "n_ZEB_SLUG" => 41,
    "fc_ZEB_ZEB" => 42,
    "fc_ZEB_SNAI" => 43,
    "fc_ZEB_SLUG" => 44,
    "K_SNAI_SNAI" => 45,
    "K_SNAI_SLUG" => 46,
    "n_SNAI_SNAI" => 47,
    "n_SNAI_SLUG" => 48,
    "fc_SNAI_SNAI" => 49,
    "fc_SNAI_SLUG" => 50,
    "K_SLUG_SLUG" => 51,
    "K_SLUG_SNAI" => 52,
    "n_SLUG_SLUG" => 53,
    "n_SLUG_SNAI" => 54,
    "fc_SLUG_SLUG" => 55,
    "fc_SLUG_SNAI" => 56,
    "K_CDH_SLUG" => 57,
    "K_CDH_SNAI" => 58,
    "K_CDH_ZEB" => 59,
    "n_CDH_SLUG" => 60,
    "n_CDH_SNAI" => 61,
    "n_CDH_ZEB" => 62,
    "fc_CDH_SLUG" => 63,
    "fc_CDH_SNAI" => 64,
    "fc_CDH_ZEB" => 65,
    "t_ZEB_2" => 66,
    "t_ZEB_22" => 67,
    "t_ZEB_1" => 68,
    "t_ZEB_12" => 69,
    "t_ZEB_122" => 70,
    "t_SLUG" => 71,
    "a_ZEB_2" => 72,
    "a_ZEB_22" => 73,
    "a_ZEB_1" => 74,
    "a_ZEB_12" => 75,
    "a_ZEB_122" => 76,
    "a_SLUG" => 77,
    "b_200_ZEB_2" => 78,
    "b_200_ZEB_22" => 79,
    "b_200_ZEB_12" => 80,
    "b_200_ZEB_122" => 81,
    "b_101_ZEB_1" => 82,
    "b_101_ZEB_12" => 83,
    "b_101_ZEB_122" => 84,
    "b_SLUG" => 85,
    "kOff" => 86,
    "KA_ZEB_101" => 87,
    "KA_ZEB_200" => 88,
    "KA_SLUG" => 89,
    "σR_101" => 90,
    "σR_200" => 91,
    "σR_ZEB" => 92,
    "σR_SNAI" => 93,
    "σR_SLUG" => 94,
    "σR_CDH" => 95,
    "σX_ZEB" => 96,
    "σX_SNAI" => 97,
    "σX_SLUG" => 98,
    "σX_CDH" => 99,
    "σC_ZEB_2" => 100,
    "σC_ZEB_22" => 101,
    "σC_ZEB_1" => 102,
    "σC_ZEB_12" => 103,
    "σC_ZEB_122" => 104,
    "σC1_SLUG" => 105,
    "τ" => 106,
)

p_default = [
    10, # 1 k_101
    10, # 2 k_200
    0.01, # 3 k_ZEB
    0.1, # 4 k_SNAI
    1, # 5 k_SLUG
    10, # 6 k_CDH
    1, # 7 d_101
    1, # 8 d_200
    1, # 9 dR_ZEB
    1, # 10 dR_SNAI
    1, # 11 dR_SLUG
    1, # 12 dR_CDH
    1, # 13 l_ZEB
    1.8, # 14 l_SNAI
    5, # 15 l_SLUG
    1, # 16 l_CDH
    1, # 17 dP_ZEB
    1.25, # 18 dP_SNAI
    1.1, # 19 dP_SLUG
    1, # 20 dP_CDH
    1.9, # 21 K_101_SNAI
    2.2, # 22 K_101_SLUG
    2, # 23 n_101_SNAI
    1, # 24 n_101_SLUG
    0.1, # 25 fc_101_SNAI
    0.4, # 26 fc_101_SLUG
    2.2, # 27 K_200_ZEB
    1.9, # 28 K_200_SNAI
    2.2, # 29 K_200_SLUG
    3, # 30 n_200_ZEB
    2, # 31 n_200_SNAI
    1, # 32 n_200_SLUG
    0.1, # 33 fc_200_ZEB
    0.1, # 34 fc_200_SNAI
    0.4, # 35 fc_200_SLUG
    0.25, # 36 K_ZEB_ZEB
    1.8, # 37 K_ZEB_SNAI
    2, # 38 K_ZEB_SLUG
    2, # 39 n_ZEB_ZEB
    2, # 40 n_ZEB_SNAI
    2, # 41 n_ZEB_SLUG
    7.5, # 42 fc_ZEB_ZEB
    10, # 43 fc_ZEB_SNAI
    4, # 44 fc_ZEB_SLUG
    1.8, # 45 K_SNAI_SNAI
    2.25, # 46 K_SNAI_SLUG
    5, # 47 n_SNAI_SNAI
    3, # 48 n_SNAI_SLUG
    0.4, # 49 fc_SNAI_SNAI
    0.5, # 50 fc_SNAI_SLUG
    2.5, # 51 K_SLUG_SLUG
    1.8, # 52 K_SLUG_SNAI
    4, # 53 n_SLUG_SLUG
    1, # 54 n_SLUG_SNAI
    4, # 55 fc_SLUG_SLUG
    0.5, # 56 fc_SLUG_SNAI
    0.5, # 57 K_CDH_SLUG
    0.5, # 58 K_CDH_SNAI
    1, # 59 K_CDH_ZEB
    2, # 60 n_CDH_SLUG
    2, # 61 n_CDH_SNAI
    2, # 62 n_CDH_ZEB
    0.3, # 63 fc_CDH_SLUG
    0.3, # 64 fc_CDH_SNAI
    0.75, # 65 fc_CDH_ZEB
    0.1, # 66 t_ZEB_2
    0.1, # 67 t_ZEB_22
    0.1, # 68 t_ZEB_1
    0.1, # 69 t_ZEB_12
    0.01, # 70 t_ZEB_122
    0.1, # 71 t_SLUG
    2, # 72 a_ZEB_2
    4, # 73 a_ZEB_22
    2, # 74 a_ZEB_1
    4, # 75 a_ZEB_12
    8, # 76 a_ZEB_122
    2, # 77 a_SLUG
    0.5, # 78 b_200_ZEB_2
    0.25, # 79 b_200_ZEB_22
    0.5, # 80 b_200_ZEB_12
    0.25, # 81 b_200_ZEB_122
    0.5, # 82 b_101_ZEB_1
    0.5, # 83 b_101_ZEB_12
    0.5, # 84 b_101_ZEB_122
    0.5, # 85 b_SLUG
    100, # 86 kOff
    10, # 87 KA_ZEB_101
    10, # 88 KA_ZEB_200
    10, # 89 KA_SLUG
    0, # 90 σR_101
    0, # 91 σR_200
    0, # 92 σR_ZEB
    0, # 93 σR_SNAI
    0, # 94 σR_SLUG
    0, # 95 σR_CDH
    0, # 96 σX_ZEB
    0, # 97 σX_SNAI
    0, # 98 σX_SLUG
    0, # 99 σX_CDH
    0, # 100 σC_ZEB_2
    0, # 101 σC_ZEB_22
    0, # 102 σC_ZEB_1
    0, # 103 σC_ZEB_12
    0, # 104 σC_ZEB_122
    0, # 105 σC1_SLUG
    0.1, # 106 τ
]
