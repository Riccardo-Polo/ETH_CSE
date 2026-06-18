from const import EstimatorConstant
import numpy as np


_NUM_STATES = 9


def _symmetrize(P: np.ndarray) -> np.ndarray:
    return 0.5 * (P + P.T)


def _theta_star(k: int, c: EstimatorConstant) -> float:
    if k < c.N * c.k_switch_frac:
        return c.theta_dive
    return -c.theta_dive


def _initial_mean_and_covariance(c: EstimatorConstant) -> tuple[np.ndarray, np.ndarray]:
    xm = np.zeros(_NUM_STATES)
    xm[2] = c.z0
    xm[4] = 0.5 * (c.theta0_min + c.theta0_max)
    xm[5] = 0.5 * (c.p0_min + c.p0_max)

    Pm = np.zeros((_NUM_STATES, _NUM_STATES))
    Pm[0, 0] = c.R ** 2 / 4.0
    Pm[1, 1] = c.R ** 2 / 4.0
    Pm[3, 3] = c.psi_bound ** 2 / 3.0
    Pm[4, 4] = (c.theta0_max - c.theta0_min) ** 2 / 12.0
    Pm[5, 5] = (c.p0_max - c.p0_min) ** 2 / 12.0
    return xm, Pm


def _process_noise_covariance(c: EstimatorConstant) -> np.ndarray:
    Q = np.zeros((_NUM_STATES, _NUM_STATES))
    Q[6, 6] = c.sigma_r_psi ** 2
    Q[7, 7] = c.sigma_r_theta ** 2
    Q[8, 8] = c.sigma_a ** 2
    return Q


def _motion_model(
        x: np.ndarray,
        k: int,
        c: EstimatorConstant,
        clip_depth: bool = False,
) -> np.ndarray:
    Ts = c.Ts
    xp = np.empty_like(x)

    psi = x[3]
    theta = x[4]
    p = x[5]
    r_psi = x[6]
    r_theta = x[7]
    a = x[8]

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)

    xp[0] = x[0] + p * cos_theta * cos_psi * Ts
    xp[1] = x[1] + p * cos_theta * sin_psi * Ts
    xp[2] = x[2] + p * sin_theta * Ts
    if clip_depth and xp[2] < 0.0:
        xp[2] = 0.0
    xp[3] = psi + r_psi * Ts
    xp[4] = theta + r_theta * Ts
    xp[5] = p + (a - c.cd * p) * Ts
    xp[6] = r_psi + c.K_yaw * (c.omega_yaw - r_psi)
    xp[7] = r_theta + c.K_theta_p * (_theta_star(k, c) - theta) - c.K_theta_d * r_theta
    xp[8] = a + c.K_p_p * (c.p_cruise - p) - c.K_p_d * a

    return xp


def _motion_jacobian(x: np.ndarray, c: EstimatorConstant) -> np.ndarray:
    Ts = c.Ts
    F = np.eye(_NUM_STATES)

    psi = x[3]
    theta = x[4]
    p = x[5]

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)

    F[0, 3] = -p * cos_theta * sin_psi * Ts
    F[0, 4] = -p * sin_theta * cos_psi * Ts
    F[0, 5] = cos_theta * cos_psi * Ts

    F[1, 3] = p * cos_theta * cos_psi * Ts
    F[1, 4] = -p * sin_theta * sin_psi * Ts
    F[1, 5] = cos_theta * sin_psi * Ts

    F[2, 4] = p * cos_theta * Ts
    F[2, 5] = sin_theta * Ts

    F[3, 6] = Ts
    F[4, 7] = Ts
    F[5, 5] = 1.0 - c.cd * Ts
    F[5, 8] = Ts
    F[6, 6] = 1.0 - c.K_yaw
    F[7, 4] = -c.K_theta_p
    F[7, 7] = 1.0 - c.K_theta_d
    F[8, 5] = -c.K_p_p
    F[8, 8] = 1.0 - c.K_p_d

    return F


def _predict(
        xm: np.ndarray,
        Pm: np.ndarray,
        k: int,
        c: EstimatorConstant,
        Q: np.ndarray,
        clip_depth: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    F = _motion_jacobian(xm, c)
    xp = _motion_model(xm, k, c, clip_depth=clip_depth)
    Pp = F @ Pm @ F.T + Q
    return xp, _symmetrize(Pp)


def _usbl_variance(depth: float, c: EstimatorConstant) -> float:
    sigma_u = c.sigma_base + c.alpha * max(0.0, depth)
    return sigma_u ** 2


def _gaussian_measurement_matrices(
        xp: np.ndarray,
        measurement: np.ndarray,
        c: EstimatorConstant,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    values = []
    variances = []

    usbl_var = _usbl_variance(xp[2], c)
    if not np.isnan(measurement[0]):
        row = np.zeros(_NUM_STATES)
        row[0] = 1.0
        rows.append(row)
        values.append(measurement[0])
        variances.append(usbl_var)
    if not np.isnan(measurement[1]):
        row = np.zeros(_NUM_STATES)
        row[1] = 1.0
        rows.append(row)
        values.append(measurement[1])
        variances.append(usbl_var)

    sonar_var = c.sigma_S ** 2
    for state_idx, meas_idx in ((0, 2), (1, 3), (2, 4)):
        if not np.isnan(measurement[meas_idx]):
            row = np.zeros(_NUM_STATES)
            row[state_idx] = 1.0
            rows.append(row)
            values.append(measurement[meas_idx])
            variances.append(sonar_var)

    if not rows:
        return (
            np.zeros((0, _NUM_STATES)),
            np.zeros(0),
            np.zeros((0, 0)),
        )

    return np.vstack(rows), np.asarray(values), np.diag(variances)


def _kalman_update(
        xp: np.ndarray,
        Pp: np.ndarray,
        measurement: np.ndarray,
        c: EstimatorConstant,
) -> tuple[np.ndarray, np.ndarray]:
    H, z, R = _gaussian_measurement_matrices(xp, measurement, c)
    if H.shape[0] == 0:
        return xp, Pp

    innovation = z - H @ xp
    S = H @ Pp @ H.T + R
    K = np.linalg.solve(S, H @ Pp).T
    xm = xp + K @ innovation
    Pm = Pp - K @ S @ K.T
    return xm, _symmetrize(Pm)


def _non_gaussian_measurement_update(
        xp: np.ndarray,
        Pp: np.ndarray,
        measurement: np.ndarray,
        c: EstimatorConstant,
) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    values = []
    variances = []
    sonar_positions = []

    usbl_var = _usbl_variance(xp[2], c)
    if not np.isnan(measurement[0]):
        row = np.zeros(_NUM_STATES)
        row[0] = 1.0
        rows.append(row)
        values.append(measurement[0])
        variances.append(usbl_var)
    if not np.isnan(measurement[1]):
        row = np.zeros(_NUM_STATES)
        row[1] = 1.0
        rows.append(row)
        values.append(measurement[1])
        variances.append(usbl_var)

    beta = c.sonar_bimodality
    component_var = max(1e-12, (1.0 - 9.0 * beta ** 2) * c.sigma_S ** 2)
    for state_idx, meas_idx in ((0, 2), (1, 3), (2, 4)):
        if not np.isnan(measurement[meas_idx]):
            row = np.zeros(_NUM_STATES)
            row[state_idx] = 1.0
            rows.append(row)
            values.append(measurement[meas_idx])
            variances.append(component_var)
            sonar_positions.append(len(rows) - 1)

    if not rows:
        return xp, Pp

    H = np.vstack(rows)
    z = np.asarray(values)
    R = np.diag(variances)
    S = H @ Pp @ H.T + R
    K = np.linalg.solve(S, H @ Pp).T
    base_cov = _symmetrize(Pp - K @ S @ K.T)

    predicted_measurement = H @ xp
    sign, logdet = np.linalg.slogdet(S)
    if sign <= 0.0:
        return _kalman_update(xp, Pp, measurement, c)

    n_meas = H.shape[0]
    n_sonar = len(sonar_positions)
    if n_sonar == 0:
        innovation = z - predicted_measurement
        return xp + K @ innovation, base_cov

    primary_mean = beta * c.sigma_S
    secondary_mean = -9.0 * beta * c.sigma_S
    primary_weight = c.sonar_mix_weight
    secondary_weight = 1.0 - primary_weight

    component_means = []
    log_weights = []
    for mask in range(1 << n_sonar):
        noise_mean = np.zeros(n_meas)
        log_prior = 0.0
        for bit, pos in enumerate(sonar_positions):
            if mask & (1 << bit):
                noise_mean[pos] = secondary_mean
                log_prior += np.log(secondary_weight)
            else:
                noise_mean[pos] = primary_mean
                log_prior += np.log(primary_weight)

        innovation = z - predicted_measurement - noise_mean
        solved = np.linalg.solve(S, innovation)
        quad = innovation @ solved
        log_likelihood = -0.5 * (n_meas * np.log(2.0 * np.pi) + logdet + quad)
        component_means.append(xp + K @ innovation)
        log_weights.append(log_prior + log_likelihood)

    component_means = np.asarray(component_means)
    log_weights = np.asarray(log_weights)
    log_weights -= np.max(log_weights)
    weights = np.exp(log_weights)
    weight_sum = np.sum(weights)
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        best = int(np.argmax(log_weights))
        return component_means[best], base_cov
    weights /= weight_sum

    xm = weights @ component_means
    Pm = base_cov.copy()
    for weight, mean in zip(weights, component_means):
        diff = mean - xm
        Pm += weight * np.outer(diff, diff)

    if xm[2] < 0.0:
        xm[2] = 0.0
    if xm[5] < 0.0:
        xm[5] = 0.0
    return xm, _symmetrize(Pm)


class EKF:
    """
    Extended Kalman Filter for the whale tracking problem under Gaussian noise.

    Args:
        estimator_constant : EstimatorConstant
            Constants known to the estimator (initial state bounds, process
            noise paramaters, etc...)
    """

    def __init__(
            self,
            estimator_constant: EstimatorConstant,
    ):
        self.constant = estimator_constant
        self.Q = _process_noise_covariance(estimator_constant)
        self.k = 0

        # Posterior mean of the state estimate
        self.xm = None
        # Posterior covariance of the state estimate
        self.Pm = None

    def initialize(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Initialize the estimator. Set up any internal state required by the
        filter and return the initial posterior mean and covariance.

        Returns:
            xm : np.ndarray, dim: (num_states,)
                The initial posterior state mean. The order of states is
                x = [x, y, z, psi, theta, p, r_psi, r_theta, a].
            Pm : np.ndarray, dim: (num_states, num_states)
                The initial posterior state covariance.
        """
        xm, Pm = _initial_mean_and_covariance(self.constant)
        self.xm = xm
        self.Pm = Pm
        self.k = 0

        return xm, Pm

    def estimate(
            self,
            measurement: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Update the estimator with a new measurement and return the posterior
        mean and covariance at the current time step.

        Args:
            measurement : np.ndarray, dim: (num_measurements,)
                Sensor measurements from time step k, z(k). The order of
                measurements is z = [z_Ux, z_Uy, z_Sx, z_Sy, z_Sz], where
                z_Ux, z_Uy are the USBL measurements (np.nan during dropout)
                and z_Sx, z_Sy, z_Sz are the sonar measurements.

        Returns:
            xm : np.ndarray, dim: (num_states,)
                The posterior state mean at time step k. The order of
                states is x = [x, y, z, psi, theta, p, r_psi, r_theta, a].
            Pm : np.ndarray, dim: (num_states, num_states)
                The posterior state covariance at time step k.
        """
        if self.xm is None or self.Pm is None:
            self.initialize()

        xp, Pp = _predict(self.xm, self.Pm, self.k, self.constant, self.Q)
        xm, Pm = _kalman_update(xp, Pp, measurement, self.constant)

        self.xm = xm
        self.Pm = Pm
        self.k += 1

        return xm, Pm


class NonGaussianEstimator:
    """
    Estimator for the whale tracking problem under non-Gaussian noise.

    You are free to choose any recursive estimation strategy you like for this
    class. The public interface only exposes a single state estimate; the
    internal representation of your filter is up to you.

    Args:
        estimator_constant : EstimatorConstant
            Constants known to the estimator.
    """

    def __init__(
            self,
            estimator_constant: EstimatorConstant,
    ):
        self.constant = estimator_constant
        self.Q = _process_noise_covariance(estimator_constant)
        self.k = 0
        self.xm = None
        self.Pm = None

    def initialize(self) -> np.ndarray:
        """
        Initialize the estimator. Set up any internal state required by the
        filter and return the initial state estimate.

        Returns:
            xm : np.ndarray, dim: (num_states,)
                The initial state estimate. The order of states is
                x = [x, y, z, psi, theta, p, r_psi, r_theta, a].
        """
        xm, Pm = _initial_mean_and_covariance(self.constant)
        self.xm = xm
        self.Pm = Pm
        self.k = 0

        return xm

    def estimate(
            self,
            measurement: np.ndarray,
    ) -> np.ndarray:
        """
        Update the estimator with a new measurement and return the posterior
        state estimate at the current time step.

        Args:
            measurement : np.ndarray, dim: (num_measurements,)
                Sensor measurements from time step k, z(k). The order of
                measurements is z = [z_Ux, z_Uy, z_Sx, z_Sy, z_Sz], where
                z_Ux, z_Uy are the USBL measurements (np.nan during dropout)
                and z_Sx, z_Sy, z_Sz are the sonar measurements.

        Returns:
            xm : np.ndarray, dim: (num_states,)
                The posterior state estimate at time step k. The order of
                states is x = [x, y, z, psi, theta, p, r_psi, r_theta, a].
        """
        if self.xm is None or self.Pm is None:
            self.initialize()

        xp, Pp = _predict(
            self.xm,
            self.Pm,
            self.k,
            self.constant,
            self.Q,
            clip_depth=True,
        )
        xm, Pm = _non_gaussian_measurement_update(
            xp,
            Pp,
            measurement,
            self.constant,
        )

        self.xm = xm
        self.Pm = Pm
        self.k += 1

        return xm
