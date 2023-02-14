import numpy as np


class LineString3D:
    """" Class for operating on 3D paths

    Akin to the 2D LineString class in the shapely package (but not nearly as
    fully featured).

    Parameters
    ----------
    coords : (N, 3) array
        3D coordinates that define the path
    """
    def __init__(self, coords):
        self.coords = coords
        self.length = self.segment_lengths().sum()

    def segment_lengths(self):
        """Lengths of each segment (defined by two adjacent points)"""
        deltas = np.diff(self.coords, axis=0)
        distances = np.sqrt((deltas ** 2).sum(axis=1))
        return distances

    def offset_of_point(self, point):
        """3D offset of a point compared to its projection to the LineString3D"""
        segments = np.zeros((self.coords.shape[0] - 1, self.coords.shape[1], 2))
        for i in range(self.coords.shape[0] - 1):
            segments[i, :, 0] = self.coords[i, :]
            segments[i, :, 1] = self.coords[i + 1, :]
        segment_vecs = np.squeeze(np.diff(segments, axis=2), axis=2)
        segment_lengths = self.segment_lengths()
        norms = segment_vecs / np.linalg.norm(segment_vecs, axis=1)[:, np.newaxis]
        starts = np.squeeze(segments[:, :, 0])
        inner_product = ((point[np.newaxis, :] - starts) * norms).sum(axis=1)
        projected = starts + norms * inner_product[:, np.newaxis]
        dist_to_proj = np.sqrt(((projected - point[np.newaxis, :]) ** 2).sum(axis=1))
        on_segment = (inner_product >= 0) & (inner_product <= segment_lengths)
        if on_segment.sum() > 0:
            # projects on to at least one line segment
            # so find which is closest, and use that projected point
            closest_on_segment = np.argmin(dist_to_proj[on_segment])
            closest_idx = np.arange(projected.shape[0])[on_segment][closest_on_segment]
            proj = projected[closest_idx]
            offset = point - proj
        else:
            # not on a line segment, so find closest point
            dist_to_line = np.sqrt(((self.coords - point[np.newaxis, :]) ** 2).sum(axis=1))
            closest_idx = np.argmin(dist_to_line)
            offset = point - self.coords[closest_idx, :]
        return offset

    def project(self, point, normalized=False):
        """Distance along this object to a point nearest the specified `point`

        Parameters
        ----------
            point : (1, 3) array
                3D coordinates of another point
            normalized : bool, default False
                Whether to normalize by this object's length

        Returns
        -------
            distance : float
                Distance along object
        """
        segments = np.zeros((self.coords.shape[0] - 1, self.coords.shape[1], 2))
        for i in range(self.coords.shape[0] - 1):
            segments[i, :, 0] = self.coords[i, :]
            segments[i, :, 1] = self.coords[i + 1, :]
        segment_vecs = np.squeeze(np.diff(segments, axis=2))
        segment_lengths = self.segment_lengths()
        norms = segment_vecs / np.linalg.norm(segment_vecs, axis=1)[:, np.newaxis]
        starts = np.squeeze(segments[:, :, 0])
        inner_product = ((point[np.newaxis, :] - starts) * norms).sum(axis=1)
        projected = starts + norms * inner_product[:, np.newaxis]
        dist_to_proj = np.sqrt(((projected - point[np.newaxis, :]) ** 2).sum(axis=1))
        on_segment = (inner_product >= 0) & (inner_product <= segment_lengths)
        if on_segment.sum() > 0:
            # projects on to at least one line segment
            # so find which is closest
            closest_on_segment = np.argmin(dist_to_proj[on_segment])
            closest_idx = np.arange(projected.shape[0])[on_segment][closest_on_segment]
            length_along_path = segment_lengths[:closest_idx].sum() + inner_product[closest_idx]
        else:
            # not on a line segment, so find closest point
            dist_to_line = np.sqrt(((self.coords - point[np.newaxis, :]) ** 2).sum(axis=1))
            closest_idx = np.argmin(dist_to_line)
            if closest_idx == 0:
                length_along_path = 0.
            else:
                length_along_path = segment_lengths[:closest_idx].sum()
        if normalized:
            return length_along_path / segment_lengths.sum()
        else:
            return length_along_path

    def rotation_to_vector(self, v):
        """Rotation matrix that will align this object to a vector `v`

        Based on https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

        Parameters
        ----------
            v : (1, 3) array
                Vector to calculate rotation toward

        Returns
        -------
            rot : (3, 3) array
                Rotation matrix
        """

        # Make a normal vector pointing from start of path to end
        path_vec = self.coords[-1, :] - self.coords[0, :]
        path_vec = path_vec / np.linalg.norm(path_vec)

        # Ensure v is unit-length
        v = v / np.linalg.norm(v)

        cross = np.cross(path_vec, v)
        dot = np.dot(path_vec, v)
        cross_mat = np.array([
            [0, -cross[2], cross[1]],
            [cross[2], 0, -cross[0]],
            [-cross[1], cross[0], 0]
        ])

        rot = np.identity(3) + cross_mat + np.dot(cross_mat, cross_mat) * (1 / (1 + dot))
        return rot
