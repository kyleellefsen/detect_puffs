"""Test for the gaussianFitting module.

python -m pytest test_gaussianFitting.py
"""
import numpy as np
import pytest
import gaussianFitting


def test_fitSingleSymmetricalGaussian():
    image_size = 30
    x = np.arange(image_size, dtype=float)
    y = np.arange(image_size, dtype=float)
    p_true = np.array([5.1, 18.5777, 2.2, 3])
    image = gaussianFitting.gaussian(x[:,None], y[None,:],
                                     xorigin=p_true[0],
                                     yorigin=p_true[1],
                                     sigma=p_true[2],
                                     amplitude=p_true[3])
    rng = np.random.RandomState(30)
    noise = 0.5 * rng.rand(image_size, image_size)
    p0 = (5, 5, 3, 1)
    bounds = [(0, 0, 1, 0), (30, 30, 10, 5)]
    p_fit, _, _ = gaussianFitting.fitGaussian(image + noise, p0, bounds, display=False)
    for i in range(4):
        print(f"p_true: {p_true[i]:.3f} ({bounds[0][i]}, {p0[i]}->{p_fit[i]:.4f}, {bounds[1][i]})")
        assert np.isclose(p_fit[i], p_true[i], atol=0.5)

def test_fitSingleSymmetricalGaussian_zeronoise():
    image_size = 30
    x = np.arange(image_size, dtype=float)
    y = np.arange(image_size, dtype=float)
    p_true = np.array([5.1, 18.5777, 2.2, 3])
    image = gaussianFitting.gaussian(x[:,None], y[None,:],
                                     xorigin=p_true[0],
                                     yorigin=p_true[1],
                                     sigma=p_true[2],
                                     amplitude=p_true[3])
    rng = np.random.RandomState(30)
    noise = 0 * rng.rand(image_size, image_size)
    p0 = (5, 5, 8, 3)
    bounds = [(0, 0, 1, 0), (30, 30, 10, 5)]
    p_fit, _, _ = gaussianFitting.fitGaussian(image + noise, p0, bounds, display=False)
    for i in range(4):
        print(f"p_true: {p_true[i]:.3f} ({bounds[0][i]}, {p0[i]}->{p_fit[i]:.4f}, {bounds[1][i]})")
        assert np.isclose(p_fit[i], p_true[i], atol=0.5)

def test_fitMultipleSymmetricalGaussian():
    image_size = 30
    x = np.arange(image_size, dtype=float)
    y = np.arange(image_size, dtype=float)
    p1_true = np.array([5.1, 18.5777, 2.2, 3])
    image1 = gaussianFitting.gaussian(x[:,None], y[None,:],
                                     xorigin=p1_true[0],
                                     yorigin=p1_true[1],
                                     sigma=p1_true[2],
                                     amplitude=p1_true[3])
    p2_true = np.array([20.2, 11.1, 3.2, 2.1])
    image2 = gaussianFitting.gaussian(x[:,None], y[None,:],
                                     xorigin=p2_true[0],
                                     yorigin=p2_true[1],
                                     sigma=p2_true[2],
                                     amplitude=p2_true[3])
    image = image1 + image2
    p_true = np.concatenate((p1_true, p2_true))
    rng = np.random.RandomState(30)
    noise = 0.5 * rng.rand(image_size, image_size)
    p0 = (5, 5, 3, 1, 8, 2, 3, 1)
    bounds = [(0, 0, 1, 0, 0, 0, 1, 0), (30, 30, 10, 5, 30, 30, 10, 5)]
    p_fit, image_fit_1, image_fit_2 = gaussianFitting.fitGaussian(image + noise, p0, bounds, nGaussians=2, display=False )
    for i in range(4):
        print(f"p_true: {p1_true[i]:.3f} ({bounds[0][i]}, {p0[i]}->{p_fit[i]:.4f}, {bounds[1][i]})")
        assert np.isclose(p_fit[i], p_true[i], atol=0.5)

def test_fitSingleRotatableGaussian_zeronoise():
    image_size = 30
    x = np.arange(image_size, dtype=float)
    y = np.arange(image_size, dtype=float)
    p_true = np.array([5.1, 18.5777, 1, 3, 45, 2.2])
    image = gaussianFitting.gaussian_rot(x[:,None], y[None,:],
                                     xorigin=p_true[0],
                                     yorigin=p_true[1],
                                     sigmax=p_true[2],
                                     sigmay=p_true[3],
                                     angle=p_true[4],
                                     amplitude=p_true[5])
    rng = np.random.RandomState(30)
    noise = 0.5 * rng.rand(image_size, image_size)
    p0 = (5, 5, 8, 8, 40, 3)
    bounds = [(0, 0, 1, 1, 0, 1), (30, 30, 10, 10, 90, 5)]
    p_fit, _, _ = gaussianFitting.fitRotGaussian(image + noise, p0, bounds, display=False)
    tolerances = [0.5, 0.5, 1, 1, 5, 0.5]
    for i in range(6):
        print(f"p_true: {p_true[i]:.3f} ({bounds[0][i]}, {p0[i]}->{p_fit[i]:.4f}, {bounds[1][i]})")
        assert np.isclose(p_fit[i], p_true[i], atol=tolerances[i])



if __name__ == '__main__':
    pytest.main()
