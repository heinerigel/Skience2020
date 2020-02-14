import matplotlib.pyplot as plt
import numpy as np
import math

import tqdm.notebook as tqdm
import weakref
from matplotlib.colors import LinearSegmentedColormap
import tqdm.notebook as tqdm
from time import sleep, time

# This for animated MCMC plot backgrounds
background_cdict = {
    "red": ((0.0, 1.0, 1.0), (1.0, 0.5, 0.5)),
    "green": ((0.0, 1.0, 1.0), (1.0, 0.5, 0.5)),
    "blue": ((0.0, 1.0, 1.0), (1.0, 0.5, 0.5)),
}
background_cmap = LinearSegmentedColormap("background_cmap", background_cdict)


# Algorithm 1: Classical Metropolis-Hastings algorithm with multivariate
# normally distributed step
#
# Steps are proposed according to N(0, I).
#
#
#
def sample_mh(target, initialModel, epsilon, ns, maxtime=None, thinning=1):

    # Assert that we are given a model vector
    assert initialModel.shape == (initialModel.size, 1)

    # Set up the variables
    m = initialModel
    X = target.misfit(m)
    accepted = 0

    # Pre-allocate samples in RAM for speed
    samples = np.empty((m.size + 1, math.ceil(ns / thinning)))

    # Create a nice progress-bar
    pbar = tqdm.trange(ns, desc="Acceptance rate: 0, progress: ")

    # If we are given a maximum time, start timer now
    if maxtime is not None:
        t_end = time() + maxtime

    # Metropolis-Hastings sampling
    try:
        for i_sample in pbar:

            # Draw a perturbation with a unit Gaussian in model space
            dm = np.random.randn(m.size, 1) * epsilon

            # Add the perturbation
            m_new = m + dm

            # Compute the new misfit
            X_new = target.misfit(m_new)

            # Evaluate acceptence criterion
            if np.exp(X - X_new) > np.random.uniform(0, 1):
                # If accepted, write new parameters and misfit over the old one
                m = m_new
                X = X_new
                # And keep track of acceptance rate
                accepted += 1

            # Every 1000 steps, update the progress bar (too often, and we
            # slow down sampling)
            if i_sample % 1000 == 0:
                pbar.set_description(
                    f"Acceptance rate: {accepted / (i_sample+1):.2}, progress: "
                )

            if thinning == 1:
                # This is important; in BOTH cases (accept or reject); write out
                # the state after the acceptance step; this is either the new or old.
                samples[:-1, i_sample] = m[:, 0]
                samples[-1, i_sample] = X
            elif i_sample % thinning == 0:
                samples[:-1, int(i_sample / thinning)] = m[:, 0]
                samples[-1, int(i_sample / thinning)] = X

            # If we are over time, interrupt
            if maxtime is not None and t_end < time():
                raise KeyboardInterrupt
    except KeyboardInterrupt:  # Catch SIGINT --------------------------------------
        # Close tqdm progressbar
        pbar.close()
        # Delete all non-written samples
        samples = samples[:, : math.ceil(i_sample / thinning) - 1]
    finally:
        # Output final acceptance rate
        pbar.set_description(
            f"Final acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
        )
        return samples


def visual_sample_mh(
    target,
    initialModel,
    epsilon,
    ns,
    extent=None,
    dims_to_visualize=[0, 1],
    fig=None,
    ax=None,
    background_samples=None,
    true_m=None,
    animate_sample_interval=10,
    animate_proposal=True,
    figsize=(6, 6),
):

    # EXTRA PLOTTING COMMANDS BELOW =============================================================================
    # Creating figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(linewidth=0)
    ax.set_xlabel(f"Model space dimension {dims_to_visualize[0]}")
    ax.set_ylabel(f"Model space dimension {dims_to_visualize[1]}")

    histogram_data, xe, ye = np.histogram2d(
        background_samples[dims_to_visualize[0], :],
        background_samples[dims_to_visualize[1], :],
        bins=int((background_samples[0, :].size) ** 0.3),
    )
    extent_histogram = [xe[0], xe[-1], ye[0], ye[-1]]
    ax.imshow(
        histogram_data.T,
        extent=extent_histogram,
        origin="lower",
        cmap=background_cmap,
    )
    if extent is None:
        extent = [[xe[0], xe[-1]], [ye[0], ye[-1]]]

    plt.xlim(extent[0])
    plt.ylim(extent[1])

    ax.scatter(
        true_m[dims_to_visualize[0]],
        true_m[dims_to_visualize[1]],
        s=40,
        color="k",
        label="True model",
    )
    ax.scatter(
        initialModel[dims_to_visualize[0]],
        initialModel[dims_to_visualize[1]],
        color="g",
        s=40,
        label="MH Samples",
    )
    plt.legend(loc=2)
    ax.set_aspect("auto")
    fig.canvas.draw()
    sleep(0.001)
    # END =======================================================================================================

    # Assert that we are given a model vector
    assert initialModel.shape == (initialModel.size, 1)

    # Set up the variables
    m = initialModel
    X = target.misfit(m)
    accepted = 0

    # Pre-allocate samples in ram for speed
    samples = np.empty((m.size + 1, ns))

    # Create a nice progress-bar
    pbar = tqdm.trange(ns, desc="Acceptance rate: 0, progress: ")

    # Metropolis-Hastings sampling
    for i_sample in pbar:

        # Draw a perturbation with a unit Gaussian in model space
        dm = np.random.randn(initialModel.size, 1) * epsilon

        # Add the perturbation
        m_new = m + dm

        # EXTRA PLOTTING COMMANDS BELOW =========================================================================
        if animate_proposal:
            proposal_line = ax.plot(
                [m[dims_to_visualize[0]], m_new[dims_to_visualize[0]]],
                [m[dims_to_visualize[1]], m_new[dims_to_visualize[1]]],
                "r",
                alpha=0.25,
            )
            plt.xlim(extent[0])
            plt.ylim(extent[1])
            fig.canvas.draw()
            fig.canvas.flush_events()
            sleep(0.01)
        # END ===================================================================================================

        # Compute the new misfit
        X_new = target.misfit(m_new)

        # Evaluate acceptence criterion
        if np.exp(X - X_new) > np.random.uniform(0, 1):
            # If accepted, write new parameters and misfit over the old one
            m = m_new
            X = X_new

            # And keep track of acceptance rate
            accepted += 1

        if i_sample % 1 == 0:
            pbar.set_description(
                f"Acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
            )

        if animate_proposal:
            l = proposal_line.pop(0)
            wl = weakref.ref(l)
            l.remove()
            del l

        # EXTRA PLOTTING COMMANDS BELOW =========================================================================
        if i_sample % animate_sample_interval == 0:
            if i_sample % 1 == 0:
                plt.cla()
                ax.set_xlabel(f"Model space dimension {dims_to_visualize[0]}")
                ax.set_ylabel(f"Model space dimension {dims_to_visualize[1]}")
                ax.imshow(
                    histogram_data.T,
                    extent=extent_histogram,
                    origin="lower",
                    cmap=background_cmap,
                )

                ax.scatter(
                    true_m[dims_to_visualize[0]],
                    true_m[dims_to_visualize[1]],
                    s=20,
                    color="k",
                    label="True model",
                )
                ax.scatter(
                    samples[dims_to_visualize[0], :i_sample],
                    samples[dims_to_visualize[1], :i_sample],
                    color="g",
                    s=40,
                    label="MH samples",
                )
                plt.legend(loc=2)

            ax.scatter(
                samples[dims_to_visualize[0], :i_sample],
                samples[dims_to_visualize[1], :i_sample],
                color="g",
                s=40,
                label="MH samples",
            )
            ax.set_aspect("auto")
            plt.xlim(extent[0])
            plt.ylim(extent[1])
            fig.canvas.draw()
            fig.canvas.flush_events()
            sleep(0.01)
        # END ===================================================================================================

        # This is important; in BOTH cases (accept or reject); write out
        # the state after the acceptance step; this is either the new or old.
        samples[:-1, i_sample] = m[:, 0]
        samples[-1, i_sample] = X

    return samples


def sample_hmc(target, initialModel, nt, dt, ns, mass_matrix):

    # Make sure the initial model is a column vector
    assert initialModel.shape == (initialModel.size, 1)

    # Create samples object. This is a lot faster than appending
    # to a list, as all the memory is requested in one go. One
    # Additional dimension is allocated for the misfits of every
    # model.
    samples = np.empty((initialModel.size + 1, ns))

    m = np.copy(initialModel)
    accepted = 0

    # Invert the mass matrix to accelerate computing. This might not always
    # be possible for large laplacian-like mass matrices. In such cases
    # sparse solves in  every leapfrog iteration are preferred.
    inv_mass = np.linalg.inv(mass_matrix)

    # Create a progress bar
    pbar = tqdm.trange(ns, desc="Acceptance rate: 0, progress: ")

    # Start sampling
    for i_sample in pbar:

        # Propose momentum
        p = np.random.multivariate_normal(
            np.zeros((m.size,)), mass_matrix, 1
        ).T

        # Randomize integration parameters a little
        nt_local = int(np.random.uniform(low=0.5, high=1.5, size=None) * nt)
        dt_local = np.random.uniform(low=0.5, high=1.5, size=None) * dt

        # Compute the initial misfit, kinetic energy and Hamiltonian
        X = target.misfit(m)
        K = 0.5 * (p.T @ (np.linalg.inv(mass_matrix) @ p)).item()

        # Implement total energy ---------------------------------------------
        H = 0
        # --------------------------------------------------------------------
        # Solution
        H = K + X

        # Create the new state as a copy
        m_new = np.copy(m)
        p_new = np.copy(p)

        # And start integrating the new state
        p_new -= (dt_local / 2.0) * target.grad(m_new)
        for i in range(nt_local - 1):
            # Implement the loop ---------------------------------------------

            # ----------------------------------------------------------------
            # Solution
            m_new += dt_local * inv_mass @ p_new
            p_new -= dt_local * target.grad(m_new)

        m_new += dt_local * inv_mass @ p_new
        p_new -= (dt_local / 2.0) * target.grad(m_new)

        # Compute the new misfit, kinetic energy and Hamiltonian
        X_new = target.misfit(m_new)
        K_new = 0.5 * (p_new.T @ inv_mass @ p_new).item()
        H_new = K_new + X_new

        # Accept if statistically relevant
        if np.exp(H - H_new) > np.random.uniform(low=0.0, high=1.0, size=None):
            m = np.copy(m_new)
            p = np.copy(p_new)
            X = X_new
            accepted += 1

        if i_sample % 100 == 0:
            pbar.set_description(
                f"Acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
            )

        # Add the new sample to the current sample. If the proposal was rejected,
        # we add the previous sample again.
        samples[:-1, i_sample] = m[:, 0]
        samples[-1, i_sample] = X

    # Output final acceptance rate
    pbar.set_description(
        f"Final acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
    )

    return samples


def visual_sample_hmc(
    target,
    initialModel,
    nt,
    dt,
    ns,
    mass_matrix,
    extent=None,
    dims_to_visualize=[0, 1],
    fig=None,
    ax=None,
    background_samples=None,
    true_m=None,
    animate_sample_interval=10,
    animate_trajectory=True,
    animate_trajectory_interval=10,
    figsize=(6, 6),
):

    # EXTRA PLOTTING COMMANDS BELOW =============================================================================
    if animate_trajectory:
        animate_sample_interval = 1
    # Creating figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(linewidth=0)
    ax.set_xlabel(f"Model space dimension {dims_to_visualize[0]}")
    ax.set_ylabel(f"Model space dimension {dims_to_visualize[1]}")

    histogram_data, xe, ye = np.histogram2d(
        background_samples[dims_to_visualize[0], :],
        background_samples[dims_to_visualize[1], :],
        bins=int((background_samples[0, :].size) ** 0.3),
    )
    extent_histogram = [xe[0], xe[-1], ye[0], ye[-1]]
    ax.imshow(
        histogram_data.T,
        extent=extent_histogram,
        origin="lower",
        cmap=background_cmap,
    )
    if extent is None:
        extent = [[xe[0], xe[-1]], [ye[0], ye[-1]]]

    plt.xlim(extent[0])
    plt.ylim(extent[1])

    ax.scatter(
        true_m[dims_to_visualize[0]],
        true_m[dims_to_visualize[1]],
        s=40,
        color="k",
        label="True model",
    )
    ax.scatter(
        initialModel[dims_to_visualize[0]],
        initialModel[dims_to_visualize[1]],
        color="g",
        s=40,
        label="HMC Samples",
    )
    plt.legend(loc=2)
    ax.set_aspect("auto")
    fig.canvas.draw()
    sleep(0.001)
    # END =======================================================================================================

    # Make sure the initial model is a column vector
    assert initialModel.shape == (initialModel.size, 1)

    # Create samples object. This is a lot faster than appending
    # to a list, as all the memory is requested in one go. One
    # Additional dimension is allocated for the misfits of every
    # model.
    samples = np.empty((initialModel.size + 1, ns))

    m = np.copy(initialModel)
    accepted = 0

    # Invert the mass matrix to accelerate computing. This might not always
    # be possible for large laplacian-like mass matrices. In such cases
    # sparse solves in  every leapfrog iteration are preferred.
    inv_mass = np.linalg.inv(mass_matrix)

    # Create a progress bar
    pbar = tqdm.trange(ns, desc="Acceptance rate: 0, progress: ")

    # Start sampling
    for i_sample in pbar:

        # Propose momentum
        p = np.random.multivariate_normal(
            np.zeros((m.size,)), mass_matrix, 1
        ).T

        # Randomize integration parameters a little
        nt_local = int(np.random.uniform(low=0.5, high=1.5, size=None) * nt)
        dt_local = np.random.uniform(low=0.5, high=1.5, size=None) * dt

        # Compute the initial misfit, kinetic energy and Hamiltonian
        X = target.misfit(m)
        K = 0.5 * (p.T @ (np.linalg.inv(mass_matrix) @ p)).item()

        # Implement total energy ---------------------------------------------
        H = 0
        # --------------------------------------------------------------------
        # Solution
        H = K + X

        # Create the new state as a copy
        m_new = np.copy(m)
        p_new = np.copy(p)

        # To draw the trajectory, we need to accumulate the points
        if animate_trajectory:
            accu = m
            accu = np.append(accu, m_new, axis=1)
            line = ax.plot(
                accu[dims_to_visualize[0], :],
                accu[dims_to_visualize[1], :],
                "r",
                alpha=0.25,
            )
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.xlim(extent[0])
            plt.ylim(extent[1])
            sleep(0.01)

        # And start integrating the new state
        p_new -= (dt_local / 2.0) * target.grad(m_new)
        for i in range(nt_local - 1):
            m_new += dt_local * inv_mass @ p_new
            p_new -= dt_local * target.grad(m_new)

            if animate_trajectory:
                accu = np.append(accu, m_new, axis=1)
            if animate_trajectory and i % animate_trajectory_interval == 0:
                l = line.pop(0)
                wl = weakref.ref(l)
                l.remove()
                del l
                line = ax.plot(
                    accu[dims_to_visualize[0], :],
                    accu[dims_to_visualize[1], :],
                    "r",
                    alpha=0.25,
                )
                fig.canvas.draw()
                fig.canvas.flush_events()
                sleep(0.01)

        if animate_trajectory:
            l = line.pop(0)
            wl = weakref.ref(l)
            l.remove()
            del l

        m_new += dt_local * inv_mass @ p_new
        p_new -= (dt_local / 2.0) * target.grad(m_new)

        # Compute the new misfit, kinetic energy and Hamiltonian
        X_new = target.misfit(m_new)
        K_new = 0.5 * (p_new.T @ inv_mass @ p_new).item()
        H_new = K_new + X_new

        # Accept if statistically relevant
        if np.exp(H - H_new) > np.random.uniform(low=0.0, high=1.0, size=None):
            m = np.copy(m_new)
            p = np.copy(p_new)
            X = X_new
            accepted += 1

        if i_sample % 1 == 0:
            pbar.set_description(
                f"Acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
            )

        # Add the new sample to the current sample. If the proposal was rejected,
        # we add the previous sample again.
        samples[:-1, i_sample] = m[:, 0]
        samples[-1, i_sample] = X

        if i_sample % animate_sample_interval == 0:

            if i_sample % 20 == 0:
                plt.cla()
                ax.imshow(
                    histogram_data.T,
                    extent=extent_histogram,
                    origin="lower",
                    cmap=background_cmap,
                )

                ax.scatter(
                    true_m[dims_to_visualize[0]],
                    true_m[dims_to_visualize[1]],
                    s=40,
                    color="k",
                    label="True model",
                )
                ax.scatter(
                    samples[dims_to_visualize[0], : i_sample + 1],
                    samples[dims_to_visualize[1], : i_sample + 1],
                    color="g",
                    s=40,
                    label="HMC Samples",
                )
                ax.set_xlabel(f"Model space dimension {dims_to_visualize[0]}")
                ax.set_ylabel(f"Model space dimension {dims_to_visualize[1]}")
                plt.legend(loc=2)
            ax.scatter(
                samples[dims_to_visualize[0], : i_sample + 1],
                samples[dims_to_visualize[1], : i_sample + 1],
                color="g",
                s=40,
            )
            ax.set_aspect("auto")
            fig.canvas.draw()
            fig.canvas.flush_events()
            sleep(0.01)

    # Output final acceptance rate
    pbar.set_description(
        f"Final acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
    )
    return samples


# Algorithm 5: Hamiltonian Monte Carlo with an optimized mass matrix.
#
#
#
#
#
def sample_hmc_opt(
    target, initialModel, nt, dt, ns, mass_matrix, maxtime=None, thinning=1
):

    # Make sure the initial model is a column vector
    assert initialModel.shape == (initialModel.size, 1)

    # Create samples object. This is a lot faster than appending
    # to a list, as all the memory is requested in one go. One
    # Additional dimension is allocated for the misfits of every
    # model.
    samples = np.empty((initialModel.size + 1, math.ceil(ns / thinning)))

    m = np.copy(initialModel)
    accepted = 0

    # Invert the mass matrix to accelerate computing. This might not always
    # be possible for large laplacian-like mass matrices. In such cases
    # sparse solves in  every leapfrog iteration are preferred.
    inv_mass = np.linalg.inv(mass_matrix)

    # Create a progress bar
    pbar = tqdm.trange(ns, desc="Acceptance rate: 0, progress: ")

    if maxtime is not None:
        t_end = time() + maxtime

    # Here lies the real power of this algorithm; optimizations according to
    # the mass matrix. If it is a diagonal, we can draw uncorrelated momentum
    # samples, which is much faster than drawing from an N-dimensional
    # correlated Gaussian.
    unit = False
    diagonal = False
    if np.count_nonzero(mass_matrix - np.diag(np.diagonal(mass_matrix))) == 0:
        if np.allclose(np.diagonal(mass_matrix), 1.0):
            unit = True
        diagonal = True
        mm_diag = np.diagonal(mass_matrix)
        mm_diag_sqrt = np.diagonal(mass_matrix)[:, np.newaxis] ** 0.5
        mm_diag_inv = 1.0 / np.diagonal(mass_matrix)[:, np.newaxis]
        assert np.all(mm_diag_inv == np.diagonal(inv_mass)[:, None])

    # Start sampling
    try:
        for i_sample in pbar:

            # Propose momentum, optimized style
            if unit:
                p = np.random.normal(0, 1, (m.size, 1))
            elif diagonal:
                p = np.random.normal(0, 1, (m.size, 1)) * mm_diag_sqrt
            else:
                p = np.random.multivariate_normal(
                    np.zeros((m.size,)), mass_matrix, 1
                ).T

            # Randomize integration parameters a little
            nt_local = int(
                np.random.uniform(low=0.5, high=1.5, size=None) * nt
            )
            dt_local = np.random.uniform(low=0.5, high=1.5, size=None) * dt

            # Compute the initial misfit, kinetic energy and Hamiltonian
            X = target.misfit(m)
            if unit:
                K = 0.5 * (p.T @ p).item()
            elif diagonal:
                K = 0.5 * (p.T @ (p * mm_diag_inv)).item()
            else:
                K = 0.5 * (p.T @ inv_mass @ p).item()
            H = K + X

            # Create the new state as a copy
            m_new = np.copy(m)
            p_new = np.copy(p)

            # And start integrating the new state

            # Integrate using a leapfrog integrator ------------------
            p_new -= (dt_local / 2.0) * target.grad(m_new)
            for i_integrate in range(nt_local - 1):
                if unit:
                    m_new += dt_local * p_new
                elif diagonal:
                    m_new += dt_local * (p_new * mm_diag_inv)
                else:
                    m_new += dt_local * inv_mass @ p_new
                p_new -= dt_local * target.grad(m_new)
            if unit:
                m_new += dt_local * p_new
            elif diagonal:
                m_new += dt_local * (p_new * mm_diag_inv)
            else:
                m_new += dt_local * inv_mass @ p_new
            p_new -= (dt_local / 2.0) * target.grad(m_new)

            # Compute the new misfit, kinetic energy and Hamiltonian
            X_new = target.misfit(m_new)
            if unit:
                K_new = 0.5 * (p_new.T @ p_new).item()
            elif diagonal:
                K_new = 0.5 * (p_new.T @ (p_new * mm_diag_inv)).item()
            else:
                K_new = 0.5 * (p_new.T @ inv_mass @ p_new).item()
            H_new = K_new + X_new

            # Accept if statistically relevant
            if np.exp(H - H_new) > np.random.uniform(
                low=0.0, high=1.0, size=None
            ):
                m = np.copy(m_new)
                p = np.copy(p_new)
                # Reusing would allow us to skip the computation of X at the beginning
                # of the loop, very useful for expensive problems. Additionally, we could
                # save the last gradient to save even more computation time, which I do
                # in my FWI work.
                X = X_new
                accepted += 1

            if i_sample % 3 == 0:
                pbar.set_description(
                    f"Acceptance rate: {accepted / (i_sample+1.0):.2f}, progress: "
                )

            if thinning == 1:
                # This is important; in BOTH cases (accept or reject); write out
                # the state after the acceptance step; this is either the new or old.
                samples[:-1, i_sample] = m[:, 0]
                samples[-1, i_sample] = X
            elif i_sample % thinning == 0:
                samples[:-1, int(i_sample / thinning)] = m[:, 0]
                samples[-1, int(i_sample / thinning)] = X

            if maxtime is not None and t_end < time():
                raise KeyboardInterrupt
    except KeyboardInterrupt:  # Catch SIGINT --------------------------------------
        # Close tqdm progressbar
        pbar.close()
        # Delete all non-written samples
        samples = samples[:, : math.ceil(i_sample / thinning) - 1]
    finally:
        return samples
    return samples
