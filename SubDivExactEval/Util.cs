using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SubDivExactEval
{
    public class Util
    {
        /// <summary>
        /// Reads the integers from the byte array.
        /// </summary>
        /// <param name="bytes"></param>
        /// <param name="readPosition"></param>
        /// <param name="integers"></param>
        public static void ReadIntegers(byte[] bytes, ref int readPosition, ref int[] integers)
        {
            for (int i = 0; i < integers.Length; i++)
            {
                integers[i] = BitConverter.ToInt32(bytes, readPosition);
                readPosition += 4;
            }
        }

        /// <summary>
        /// Reads a single integer from teh byte array
        /// </summary>
        /// <param name="bytes"></param>
        /// <param name="readPosition"></param>
        /// <returns></returns>
        public static int ReadInteger(byte[] bytes, ref int readPosition)
        {
            int result = BitConverter.ToInt32(bytes, readPosition);
            readPosition += 4;
            return result;
        }

        /// <summary>
        /// Reads an array of doubles from a byte array from the given start position.
        /// </summary>
        /// <param name="bytes"></param>
        /// <param name="readPosition"></param>
        /// <param name="doubles"></param>
        public static void ReadDoubles(byte[] bytes, ref int readPosition, ref double[] doubles)
        {
            for (int i = 0; i < doubles.Length; i++)
            {
                doubles[i] = BitConverter.ToDouble(bytes, readPosition);
                readPosition += 8;
            }
        }

        /// <summary>
        /// Reads a single double from a byte array
        /// </summary>
        /// <param name="bytes"></param>
        /// <param name="readPosition"></param>
        /// <returns></returns>
        public static double ReadDouble(byte[] bytes, ref int readPosition)
        {
            double result = BitConverter.ToDouble(bytes, readPosition);
            readPosition += 8;
            return result;
        }
    }
}
