
import numpy as np
import math

def regular_grid(limits=(0, 10, 0, 10), dx=1.0, dy=1.0, sigma=0.1):
    
    # Unpack limits
    xmin, xmax, ymin, ymax = limits
    nx = 1 + (xmax - xmin) / float(dx)
    ny = 1 + (ymax - ymin) / float(dy)
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(x, y)
    points = np.asarray((xx.flatten(), yy.flatten())).T
    noise = np.random.multivariate_normal(mean=(0, 0), cov=sigma*np.eye(2), size=int(nx*ny))
    points = points + noise
    return points

    
def triangular_grid(limits=(0, 10, 0, 10), side=1, sigma=0.1):
    # Unpack limits
    xmin, xmax, ymin, ymax = limits
    nx = 1 + (xmax - xmin) / float(side)
    h = side * math.sin(math.pi / 3)
    ny = 1 + (ymax - ymin) / h
    # Generate array first
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(x, y)
    # Now we need to shift the x coordinate of every second row
    xx[1::2,] += 1/2. * side
    points = np.asarray((xx.flatten(), yy.flatten())).T
    noise = np.random.multivariate_normal(mean=(0, 0), cov=sigma*np.eye(2), size=(points.shape[0]))
    points = points + noise
    return points


def hexagonal_grid(limits=(0, 10, 0, 10), side=1, sigma=0.1):
    # Unpack limits
    xmin, xmax, ymin, ymax = limits
    nx = 1 + (xmax - xmin) / float(side)
    h = side * math.sin(math.pi / 3)
    ny = 1 + (ymax - ymin) / h
    # Generate array first
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(x, y)
    # Now we need to shift the x coordinate of every second row
    xx[1::2,] += 1/2. * side
    # This produces a triangular grid. We now remove the entries corresponding to the 
    # centers of the hexagons to obtain a hexagonal grid.
    xx[0::2, 1::3] = -1
    yy[0::2, 1::3] = -1
    xx[1::2, 2::3] = -1
    yy[1::2, 2::3] = -1
    xx = xx[xx >= 0]
    yy = yy[yy >= 0]
    points = np.asarray((xx.flatten(), yy.flatten())).T
    noise = np.random.multivariate_normal(mean=(0, 0), cov=sigma*np.eye(2), size=(points.shape[0]))
    points = points + noise
    return points


def sample_ppp(limits=(0, 1, 0, 1), rate=50):
    """
    Sample from Poisson point process in unit square (will be rescaled afterwards).
    """
    xmin, xmax, ymin, ymax = limits
    N = np.random.poisson(lam=rate)
    x = np.random.uniform(xmin, xmax, N)
    y = np.random.uniform(ymin, ymax, N)
    points = np.array([x, y]).T
    
    return points


def sample_thomas(limits=(0, 1, 0, 1), kappa=5, mu=10, sigma=0.05):
    # Thomas point process
    # sigma is the variance for the inner Gaussian
    parent_points = sample_ppp(limits=limits, rate=kappa)
    points = []
    # Sample a bivariate normal around each parent point.
    for parent in parent_points:
        # Number of child points
        m = np.random.poisson(lam=mu)
        child_points = np.random.multivariate_normal(mean=parent, cov=sigma*np.eye(2), size=m)
        points.extend(child_points)
    return np.array(points)


def sample_matern(limits=(0, 1, 0, 1), kappa=5, mu=10, r=0.3):
    """
    Sample from a Matern point process.
    """
    # First generate the parent points.
    parent_points = sample_ppp(limits=limits, rate=kappa)
    points = []
    # For each of these parent points, we now sample m ~ Poisson(mu) points uniformly in a circle of 
    # radius r around it.
    for parent in parent_points:
        # Number of child points
        m = np.random.poisson(lam=mu)
        # Radii
        radii = np.random.uniform(0, r, m)
        # Thetas
        thetas = np.random.uniform(0, 2*math.pi, m)
        # Convert to cartesian
        children = np.array([radii * np.cos(thetas), radii * np.sin(thetas)]).T + parent
        points.extend(children)
    return np.array(points)

def periodic_bc(points, limit=10.0):
    points_restricted = np.empty(points.shape)
    points_restricted = np.where(points < 0, points + limit, points)
    points_restricted = np.where(points_restricted > limit,
                                 points_restricted - limit, points_restricted)
    return points_restricted

def get_flat_torus_dists(x1, y1, x2, y2):
    """
    Compute all pairwise distances between all points (x1, y1) and points (x2, y2)
    on the flat torus [0, 1] x [0, 1]
    
    Parameters:
    x1 : ndarray (M)
        An M-length list of x coordinates of each point in the first point set
    y1 : ndarray (M)
        An M-length list of y coordinates of each point in the first point set
    x2 : ndarray (N)
        An N-length list of x coordinates of each point in the second point set
    y2 : ndarray (N)
        An N-length list of y coordinates of each point in the second point set

    Returns:
    D : ndarray (M, N)
        A distance matrix whose ijth entry is the distance along the flat torus between (x1[i], y1[i]) and (x2[j], y2[j])
    """
    
    dx = np.minimum.reduce([np.abs(x1[:, None] - x2[None, :]), 
                    np.abs((x1[:, None] + 1) - x2[None, :]),
                    np.abs(x1[:, None] - (x2[None, :] + 1))])
    dy = np.minimum.reduce([np.abs(y1[:, None] - y2[None, :]),
                    np.abs((y1[:, None] + 1) - y2[None, :]),
                    np.abs(y1[:, None] - (y2[None, :] + 1))])

    return np.sqrt(dx**2 + dy**2)