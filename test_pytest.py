import pytest
from main import Deffect, Problem, Solve

PRECISION = 1e-3

@pytest.fixture
def deffect():
    t = 34e-3
    a = 0.2 * t
    c = 0.5 * t
    return Deffect(a, c)

@pytest.fixture
def problem(deffect):
    # data for problem object
    Dout = 828e-3
    t = 34e-3
    p = 7
    T = 285
    M = 4e6
    Nnue = 300

    C = 1.5e-10
    m = 3.1
    return Problem(deffect, t, Dout, p, T, M, C, m)

@pytest.fixture
def solve_ring(deffect, problem):
    return Solve(problem, deffect, "Кольцевой эффект")

@pytest.fixture
def solve_long(deffect, problem):
    return Solve(problem, deffect, "Продольный эффект")

def test_DEFFECT_double_c(deffect, problem):
    t = 34e-3
    a = 6.8e-3
    c = 0.034 * 0.5
    assert deffect.double_c == c * 2

def test_PROBLEM_Rout(deffect, problem):
    assert problem.Rout == pytest.approx(0.414, PRECISION)

def test_PROBLEM_Rin(deffect, problem):
    assert problem.Rin == pytest.approx(0.380, PRECISION)

def test_PROBLEM_ro(deffect, problem):
    assert problem.ro == 7775

def test_SOLVE_gamma_d(solve_ring, solve_long):
    assert solve_ring.gamma_d() == pytest.approx(0.705, PRECISION)
    assert solve_long.gamma_d() == pytest.approx(0.705, PRECISION)

def test_SOLVE_m_c(solve_ring, solve_long):
    assert solve_ring.m_c() == pytest.approx(0.978, PRECISION)
    assert solve_long.m_c() == pytest.approx(0.966, PRECISION)

def test_SOLVE_m_d(solve_ring, solve_long):
    assert solve_ring.m_d() == pytest.approx(1, PRECISION)
    assert solve_long.m_d() == pytest.approx(1, PRECISION)

def test_SOLVE_calc_Y(solve_ring, solve_long):
    assert solve_ring.calc_Y('C') == pytest.approx(1.708, PRECISION)
    assert solve_ring.calc_Y('D') == pytest.approx(1.230, PRECISION)
    assert solve_long.calc_Y('C') == pytest.approx(1.686, PRECISION)
    assert solve_long.calc_Y('D') == pytest.approx(1.230, PRECISION)

def test_SOLVE_sig_m(solve_ring, solve_long):
    assert solve_ring.sig_m() == pytest.approx(37.443, PRECISION)
    assert solve_long.sig_m() == pytest.approx(78.235, PRECISION)

def test_SOLVE_sig_b(solve_ring, solve_long):
    assert solve_ring.sig_b() == pytest.approx(8.078, PRECISION)
    assert solve_long.sig_b() == pytest.approx(0)

def test_SOLVE_sig_a(solve_ring, solve_long):
    assert solve_ring.sig_a() == pytest.approx(44.858, PRECISION)
    assert solve_long.sig_a() == pytest.approx(78.235, PRECISION)

def test_SOLVE_sig_c(solve_ring, solve_long):
    assert solve_ring.sig_c() == pytest.approx(44.99, PRECISION)
    assert solve_long.sig_c() == pytest.approx(78.235, PRECISION)

def test_SOVLE_calc_Sig_eq(solve_ring, solve_long):
    assert solve_ring.calc_Sig_eq('C') == pytest.approx(44.935, PRECISION)
    assert solve_ring.calc_Sig_eq('D') == pytest.approx(44.882, PRECISION)
    assert solve_long.calc_Sig_eq('C') == pytest.approx(78.235, PRECISION)
    assert solve_long.calc_Sig_eq('D') == pytest.approx(78.235, PRECISION)

def test_SOLVE_calc_K(solve_ring, solve_long, deffect, problem):
    assert solve_ring.calc_K(solve_ring.calc_Y('C'),
                             solve_ring.calc_Sig_eq('C'),
                             deffect.a) == pytest.approx(6.328, PRECISION)
    assert solve_ring.calc_K(solve_ring.calc_Y('D'),
                             solve_ring.calc_Sig_eq('D'),
                             deffect.a) == pytest.approx(4.553, PRECISION)
    assert solve_long.calc_K(solve_long.calc_Y('C'),
                             solve_long.calc_Sig_eq('C'),
                             deffect.a) == pytest.approx(10.879, PRECISION)
    assert solve_long.calc_K(solve_long.calc_Y('D'),
                             solve_long.calc_Sig_eq('D'),
                             deffect.a) == pytest.approx(7.937, PRECISION)

def test_SOVLE_solvation(solve_ring, solve_long):
    assert solve_ring.changing_per_iter() == pytest.approx({'C': 4.573e-8,
                                                            'D': 1.648e-8,
                                                            'delK_C': 6.329,
                                                            'delK_D': 4.553}, PRECISION)
    assert solve_long.changing_per_iter() == pytest.approx({'C': 2.452e-7,
                                                            'D': 9.228e-8,
                                                            'delK_C': 10.879,
                                                            'delK_D': 7.937}, PRECISION)
