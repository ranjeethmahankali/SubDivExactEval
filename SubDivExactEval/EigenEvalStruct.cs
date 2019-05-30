using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

namespace SubDivExactEval
{
    /// <summary>
    /// Different properties of the limit surface patch that can be evaluated.
    /// </summary>
    public enum EvaluationType
    {
        /// <summary>
        /// To evaluate the position of point at the given u, v parameters.
        /// </summary>
        Position = 1,
        /// <summary>
        /// Evaluate the normal limit surface at the given u, v parameters.
        /// </summary>
        Normal = 2,
        /// <summary>
        /// Evaluates the curvature of the limit surface at the given u, v parameters.
        /// </summary>
        Curvature = 3,
        /// <summary>
        /// Evaluates the tangent in the u direction, at the given u, v parameters.
        /// </summary>
        UTangent = 4,
        /// <summary>
        /// Evaluates the tangent in the v direction, at the given u, v parameters.
        /// </summary>
        VTangent = 5
    }
    /// <summary>
    /// Class to read the eigen structures from a file that stores the precomputed values.
    /// </summary>
    public class EigenEvalStructure
    {
        public const double MINIMUM_UV = 1e-12;
        /// <summary>
        /// This is the binary file containing the precomputed eigen structures.
        /// </summary>
        private static readonly string fileName = "ccdata50NT.dat";
        private static Dictionary<int, EigenEvalStructure> _precomputedEigenStructures = new Dictionary<int, EigenEvalStructure>();
        /// <summary>
        /// Reads the eigen structures from the given code.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        static EigenEvalStructure()
        {
            string fullPath = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), fileName);
            byte[] bytes = File.ReadAllBytes(fullPath);
            int readPosition = 0;
            int nEval = Util.ReadInteger(bytes, ref readPosition);
            int n, k;
            double[] vals, vecs;
            double[][] coeffs;
            for (int i = 0; i < nEval - 2; i++)
            {
                n = i + 3;
                k = 2 * n + 8;
                vals = new double[k];
                Util.ReadDoubles(bytes, ref readPosition, ref vals);
                vecs = new double[k * k];
                Util.ReadDoubles(bytes, ref readPosition, ref vecs);
                coeffs = new double[3][];
                for (int kk = 0; kk < 3; kk++)
                {
                    coeffs[kk] = new double[k * 16];
                    Util.ReadDoubles(bytes, ref readPosition, ref coeffs[kk]);
                }

                var precomputed = new EigenEvalStructure(vals, vecs, coeffs);
                _precomputedEigenStructures[precomputed._valence] = precomputed;
            }
        }

        #region fields and properties
        private int _numControlPoints, _valence;
        private double[] _eigenValues;
        private double[] _inverseEigenVectorMatrix;
        private double[][] _coefficients;
        #endregion

        #region constructors
        private EigenEvalStructure(double[] eigenValues, double[] inverseEigenVectors, double[][] bicubicCoefficients)
        {
            _numControlPoints = eigenValues.Length;
            _valence = (_numControlPoints - 8) / 2;
            // Validate all the array sizes first;
            if (inverseEigenVectors.Length != _numControlPoints * _numControlPoints)
            {
                throw new InvalidOperationException("The eigen vector data is incomplete.");
            }

            if (bicubicCoefficients.Length != 3)
            {
                throw new InvalidOperationException("The size of the coefficient array is incorrect.");
            }

            for (int k = 0; k < 3; k++)
            {
                if (bicubicCoefficients[k].Length != _numControlPoints * 16)
                {
                    throw new InvalidOperationException("The number of coefficients is incorrect.");
                }
            }

            // Copy the eigen values and eigen vectors.
            _eigenValues = eigenValues;
            _inverseEigenVectorMatrix = inverseEigenVectors;

            // Copy the bicubic spline coefficients.
            _coefficients = new double[3][];
            for (int kki = 0; kki < _coefficients.Length; kki++)
            {
                _coefficients[kki] = new double[_numControlPoints * 16];
                for (int pi = 0; pi < _numControlPoints * 16; pi++)
                {
                    _coefficients[kki][pi] = bicubicCoefficients[kki][pi];
                }
            }
        }
        #endregion

        #region methods
        /// <summary>
        /// Checks to make sure that all the precomputed eigen structures for valences in the range [3, 50] are correctly loaded.
        /// </summary>
        /// <returns></returns>
        public static bool LoadedSuccessfully()
        {
            for (int valence = 3; valence < 51; valence++)
            {
                if (!_precomputedEigenStructures.ContainsKey(valence))
                {
                    return false;
                }
            }
            return true;
        }

        private static int FlatIndex(int i, int j, int n)
        {
            return i + (j * n);
        }

        private double[] ProjectPointsToEigenSpace(double[] controlPointCoordinates)
        {
            if (controlPointCoordinates.Length != 3 * _numControlPoints)
            {
                throw new InvalidDataException("Incorrect number of control points for the eigen structure.");
            }
            double[] projectedCoords = new double[controlPointCoordinates.Length];
            double sum;
            for (int c = 0; c < 3; c++)
            {
                for (int i = 0; i < _numControlPoints; i++)
                {
                    sum = 0;
                    for (int j = 0; j < _numControlPoints; j++)
                    {
                        sum += _inverseEigenVectorMatrix[FlatIndex(i, j, _numControlPoints)] * controlPointCoordinates[FlatIndex(j, c, _numControlPoints)];
                    }
                    projectedCoords[FlatIndex(i, c, _numControlPoints)] = sum;
                }
            }

            return projectedCoords;
        }

        private static double[] EvalSpline(double[] uTerms, double[] vTerms)
        {
            double[] bSpline = new double[16];
            for (int v = 0; v < 4; v++)
            {
                for (int u = 0; u < 4; u++)
                {
                    bSpline[v * 4 + u] = uTerms[u] * vTerms[v];
                }
            }
            return bSpline;
        }

        /// <summary>
        /// To be used in evaluating a point. P = f(u,v)
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSpline(double u, double v)
        {
            double u2, u3, v2, v3;
            u2 = u * u;
            u3 = u2 * u;
            v2 = v * v;
            v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (1 - 3 * u + 3 * u2 - u3) / 6,
                    (4 - 6 * u2 + 3 * u3) / 6,
                    (1 + 3 * u + 3 * u2 - 3 * u3) / 6,
                    u3 / 6,
                },
                new double[] {
                    (1 - 3 * v + 3 * v2 - v3) / 6,
                    (4 - 6 * v2 + 3 * v3) / 6,
                    (1 + 3 * v + 3 * v2 - 3 * v3) / 6,
                    v3 / 6
                }
            );
        }

        /// <summary>
        /// To be used in evaluating the tangents in u direction i.e. (d/du) of the surface function f(u,v).
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSplineU(double u, double v)
        {
            double u2, u3, v2, v3;
            u2 = u * u;
            u3 = u2 * u;
            v2 = v * v;
            v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (-3 + 6 * u - 3 * u2) / 6,
                    (-12 * u + 9 * u2) / 6,
                    (3 + 6 * u - 9 * u2) / 6,
                    u2 / 2
                },
                new double[] {
                    (1 - 3 * v + 3 * v2 - v3) / 6,
                    (4 - 6 * v2 + 3 * v3) / 6,
                    (1 + 3 * v + 3 * v2 - 3 * v3) / 6,
                    v3 / 6
                }
            );
        }

        /// <summary>
        /// To be used in evaluating the tangents in v direction i.e. (d/dv) of the surface function f(u,v).
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSplineV(double u, double v)
        {
            double u2, u3, v2, v3;
            u2 = u * u;
            u3 = u2 * u;
            v2 = v * v;
            v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (1 - 3 * u + 3 * u2 - u3) / 6,
                    (4 - 6 * u2 + 3 * u3) / 6,
                    (1 + 3 * u + 3 * u2 - 3 * u3) / 6,
                    u3 / 6
                },
                new double[] {
                    (-3 + 6 * v - 3 * v2) / 6,
                    (-12 * v + 9 * v2) / 6,
                    (3 + 6 * v - 9 * v2) / 6,
                    v2 / 2
                }
            );
        }

        /// <summary>
        /// To be used in computing the double derivative (curvature) d/du^2, i.e. d/du of d/du of f(u,v).
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSplineUU(double u, double v)
        {
            double u2, u3, v2, v3;
            u2 = u * u;
            u3 = u2 * u;
            v2 = v * v;
            v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (6 - 6 * u) / 6,
                    (-12 + 18 * u) / 6,
                    (6 - 18 * u) / 6,
                    u
                },
                new double[] {
                    (1 - 3 * v + 3 * v2 - v3) / 6,
                    (4 - 6 * v2 + 3 * v3) / 6,
                    (1 + 3 * v + 3 * v2 - 3 * v3) / 6,
                    v3 / 6
                }
            );
        }

        /// <summary>
        /// To be used in computing the double derivative (curvature) d/dudv, i.e. d/du of d/dv of f(u,v).
        /// Note that this double derivative is the same if you change the order of partial derivatives. i.e. d/dudv is the same as d/dvdu.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSplineUV(double u, double v)
        {
            double u2, u3, v2, v3;
            u2 = u * u;
            u3 = u2 * u;
            v2 = v * v;
            v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (-3 + 6 * u - 3 * u2) / 6,
                    (-12 * u + 9 * u2) / 6,
                    (3 + 6 * u - 9 * u2) / 6,
                    u2 / 2
                },
                new double[] {
                    (-3 + 6 * v - 3 * v2) / 6,
                    (-12 * v + 9 * v2) / 6,
                    (3 + 6 * v - 9 * v2) / 6,
                    v2 / 2
                }
            );
        }

        /// <summary>
        /// To be used in computing the double derivative (curvature) d/dv^2, i.e. d/dv of d/dv of f(u,v).
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] EvaluateBSplineVV(double u, double v)
        {
            double u2, u3, v2, v3;
            double N0u, N1u, N2u, N3u, N0v, N1v, N2v, N3v;
            double[] coeff = new double[16];

            u2 = u * u; u3 = u2 * u;
            v2 = v * v; v3 = v2 * v;

            return EvalSpline(
                new double[] {
                    (1 - 3 * u + 3 * u2 - u3) / 6,
                    (4 - 6 * u2 + 3 * u3) / 6,
                    (1 + 3 * u + 3 * u2 - 3 * u3) / 6,
                    u3 / 6
                },
                new double[] {
                    (6 - 6 * v) / 6,
                    (-12 + 18 * v) / 6,
                    (6 - 18 * v) / 6,
                    v
                }
            );
        }

        /// <summary>
        /// This is the method used for evaluating surfaces. This should not be used from anywhere else. Please use the other wrapper methods instead.
        /// </summary>
        /// <param name="controlPointCoords">Array of 2*N+8 control points where N is the valence of
        /// this patch, sorted in the order described in the Stam's paper. The
        /// coordinates of the control points should be packed into the flat array
        /// in the form - [x1, x2, x3, x4, ...y1, y2, y3, y4, ...z1, z2, z3, z4...]</param>
        /// <param name="uvParams">This should be a nested array, representing the u,v coordinates of the
        /// points to be evalutated.The inner arrays are expected to be of length 2.</param>
        /// <param name="evalType">This enumeration tells the function what property to evaluate at the
        /// given parameters.See the documentation of this enumeration for more details.</param>
        /// <param name="uvTolerance">This will be used as the minimum value of u and v since the evaluation
        /// at an extraordinary points is not possible(see the paper). If not
        /// supplied, then the default value is used defined as a const member
        /// variable(1e-12).</param>
        /// <returns></returns>
        public static double[] EvaluateSurface(double[] controlPointCoords, double[][] uvParams, EvaluationType evalType, double uvTolerance = MINIMUM_UV)
        {
            if (controlPointCoords.Length % 3 != 0)
            {
                throw new ArgumentException("The controlpoint data is invalid.");
            }
            int numControlPoints = controlPointCoords.Length / 3;
            if ((numControlPoints - 8) % 2 != 0)
            {
                throw new ArgumentException("The number of control points doesn't make sense");
            }
            int valence = (numControlPoints - 8) / 2;
            if (valence < 3 || valence > 50)
            {
                throw new ArgumentOutOfRangeException("The supplied valence value is not supported. The valence of the patch must be between " +
                    "3 and 50 (inclusive of both)");
            }

            // Sanitize the inputs first.
            if (uvTolerance < 0)
            {
                uvTolerance = -uvTolerance;
            }
            else if (uvTolerance == 0)
            {
                uvTolerance = MINIMUM_UV;
            }

            int numDerivative = 1;
            switch (evalType)
            {
                case EvaluationType.Curvature:
                    numDerivative = 3;
                    break;
                case EvaluationType.UTangent:
                case EvaluationType.VTangent:
                case EvaluationType.Normal:
                    numDerivative = 2;
                    break;
                default:
                case EvaluationType.Position:
                    numDerivative = 1;
                    break;
            }

            EigenEvalStructure eigen = _precomputedEigenStructures[valence];
            double u, v, pow2, logU, logV, u0, v0, f, splineTerm, powL;
            int coeffIndex = -1;
            double[][] splineCoefficients = new double[3][];
            // This array will either contain the coordinates of evaluated points or evaluated normals, depending on the what the input parameter asks for.
            double[] evalResults = numDerivative == 3 ? new double[uvParams.Length * 9] : new double[uvParams.Length * 3];
            double[] projectedCoordinates = eigen.ProjectPointsToEigenSpace(controlPointCoords);
            int subLevel;
            double[] tempVec;
            for (int uvi = 0; uvi < uvParams.Length; uvi++)
            {
                u = uvParams[uvi][0];
                v = uvParams[uvi][1];
                // u and v cannot be zero. They can be at the least the tolerance.
                u = Math.Max(u, uvTolerance);
                v = Math.Max(v, uvTolerance);

                logU = -Math.Log(u, 2) + 1;
                logV = -Math.Log(v, 2) + 1;

                // Determining which level of subdivision is needed for the given uv point to be inside a local regular patch.
                // This sub level would be infinity if the uv is (0,0) - coincident with the extraordinary point, which is why we ensured earlier they were not zeros.
                subLevel = logU <= logV ? (int)logU : (int)logV;

                pow2 = Math.Pow(2, subLevel - 1);
                u0 = u * pow2;
                v0 = v * pow2;

                int k;
                if (v0 < 0.5)
                {
                    coeffIndex = 0;
                    u0 = 2 * u0 - 1;
                    v0 = 2 * v0;
                }
                else if (u0 < 0.5)
                {
                    coeffIndex = 2;
                    u0 = 2 * u0;
                    v0 = 2 * v0 - 1;
                }
                else
                {
                    coeffIndex = 1;
                    u0 = 2 * u0 - 1;
                    v0 = 2 * v0 - 1;
                }

                f = 0;
                switch (numDerivative)
                {
                    case 1:
                        splineCoefficients[0] = EvaluateBSpline(u0, v0);
                        f = 1;
                        break;
                    case 2:
                        splineCoefficients[0] = EvaluateBSplineU(u0, v0);
                        splineCoefficients[1] = EvaluateBSplineV(u0, v0);
                        f = 2;
                        break;
                    case 3:
                        splineCoefficients[0] = EvaluateBSplineUU(u0, v0);
                        splineCoefficients[1] = EvaluateBSplineUV(u0, v0);
                        splineCoefficients[2] = EvaluateBSplineVV(u0, v0);
                        f = 4;
                        break;
                    default:
                        splineCoefficients[0] = EvaluateBSpline(u0, v0);
                        f = 1;
                        break;
                }

                double[][] evalData = new double[3][];
                for (int i = 0; i < 3; i++)
                {
                    evalData[i] = new double[] { 0, 0, 0 };
                }

                // Now evaluating the surface with the fancy formula.
                for (int i = 0; i < eigen._numControlPoints; i++)
                {
                    for (int d = 0; d < numDerivative; d++)
                    {
                        splineTerm = 0;
                        for (int j = 0; j < 16; j++)
                        {
                            splineTerm += eigen._coefficients[coeffIndex][FlatIndex(i, j, eigen._numControlPoints)] * splineCoefficients[d][j];
                        }
                        powL = subLevel == 0 ? 1 : subLevel == 1 ? f : f * Math.Pow(f * eigen._eigenValues[i], (double)subLevel - 1);

                        for (int c = 0; c < 3; c++)
                        {
                            evalData[d][c] += powL * projectedCoordinates[FlatIndex(i, c, eigen._numControlPoints)] * splineTerm;
                        }
                    }
                }

                // Now packing the evaluated data in a way that makes sense.
                if (numDerivative == 1)
                {
                    for (int evi = 0; evi < 3; evi++)
                    {
                        evalResults[uvi * 3 + evi] = evalData[0][evi];
                    }
                }
                else if (numDerivative == 2)
                {
                    if (evalType == EvaluationType.UTangent)
                    {
                        for (int evi = 0; evi < 3; evi++)
                        {
                            evalResults[uvi * 3 + evi] = evalData[0][evi];
                        }
                    }
                    else if (evalType == EvaluationType.VTangent)
                    {
                        for (int evi = 0; evi < 3; evi++)
                        {
                            evalResults[uvi * 3 + evi] = evalData[1][evi];
                        }
                    }
                    else
                    {
                        // The algo. gives us the tangent along u and v directions.
                        // And with this cross product we get at the normal of the surface.
                        tempVec = VectorCrossProduct(evalData[0], evalData[1]);
                        for (int cci = 0; cci < 3; cci++)
                        {
                            evalResults[uvi * 3 + cci] = tempVec[cci];
                        }
                    }
                }
                else if (numDerivative == 3)
                {
                    tempVec = evalData.SelectMany(ev => ev).ToArray();
                    for (int cci = 0; cci < 9; cci++)
                    {
                        evalResults[uvi * 9 + cci] = tempVec[cci];
                    }
                }
            }

            //return eigen.ProjectPointsToRealSpace(evalPoints);
            return evalResults;
        }

        private static double[] VectorCrossProduct(double[] a, double[] b)
        {
            return new double[] {
                a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]
            };
        }
        #endregion
    }
}
