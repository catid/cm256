#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4530) // warning C4530: C++ exception handler used, but unwind semantics are not enabled. Specify /EHsc
#endif

#include <iostream>
using namespace std;

#ifdef _MSC_VER
    #pragma warning(pop)
#endif


#include "../cm256.h"

#include <Windows.h>

static double GetPerfFrequencyInverse()
{
    // Static value containing frequency inverse of performance counter
    static double PerfFrequencyInverse = 0.;

    // If not initialized,
    if (PerfFrequencyInverse == 0.)
    {
        // Initialize the inverse (same as in Timer code)
        LARGE_INTEGER freq;
        ::QueryPerformanceFrequency(&freq);
        PerfFrequencyInverse = 1.0 / (double)freq.QuadPart;
    }

    return PerfFrequencyInverse;
}

void initializeBlocks(unsigned char* originals[256], int blockCount, int blockBytes)
{
    for (int i = 0; i < blockCount; ++i)
    {
        for (int j = 0; j < blockBytes; ++j)
        {
            const uint8_t expected = (uint8_t)(i + j * 13);
            originals[i][j] = expected;
        }
    }
}

bool validateSolution(cm256_block_t* blocks, int blockCount, int blockBytes)
{
    uint8_t seen[256] = { 0 };

    for (int i = 0; i < blockCount; ++i)
    {
        uint8_t index = blocks[i].Index;

        if (index >= blockCount)
        {
            return false;
        }

        if (seen[index])
        {
            return false;
        }

        seen[index] = 1;

        for (int j = 0; j < blockBytes; ++j)
        {
            const uint8_t expected = (uint8_t)(index + j * 13);
            if (blocks[i].Block[j] != expected)
            {
                return false;
            }
        }
    }

    return true;
}


int main()
{
    if (cm256_init())
    {
        return 1;
    }

    static const int MaxBlockBytes = 10000; // multiple of 10

    unsigned char* originals[256];

    unsigned char* orig_data = new unsigned char[256 * MaxBlockBytes];
    for (int i = 0; i < 256; ++i)
    {
        originals[i] = orig_data + i * MaxBlockBytes;
    }

    unsigned char* recoveryData = new unsigned char[256 * MaxBlockBytes];

    cm256_block_t blocks[256];

    for (int blockBytes = 8 * 162; blockBytes <= MaxBlockBytes; blockBytes *= 10)
    {
        for (int originalCount = 1; originalCount < 256; ++originalCount)
        {
            for (int recoveryCount = 1; recoveryCount <= 1 + originalCount / 2 && recoveryCount <= 256 - originalCount; ++recoveryCount)
            {
                cm256_encoder_params params;
                params.BlockBytes = blockBytes;
                params.OriginalCount = originalCount;
                params.RecoveryCount = recoveryCount;

                initializeBlocks(originals, originalCount, blockBytes);

                {
                    LARGE_INTEGER t0; ::QueryPerformanceCounter(&t0);

                    if (cm256_encode(params, originals, recoveryData))
                    {
                        cout << "Encoder error" << endl;
                        return 1;
                    }

                    LARGE_INTEGER t1; ::QueryPerformanceCounter(&t1);

                    LARGE_INTEGER tsum;
                    tsum.QuadPart = t1.QuadPart - t0.QuadPart;
                    double opusec = tsum.QuadPart * GetPerfFrequencyInverse() * 1000000.;
                    double mbps = (params.BlockBytes * params.OriginalCount / opusec);

                    cout << "Encoder: " << blockBytes << " bytes k = " << originalCount << " m = " << recoveryCount << " : " << opusec << " usec, " << mbps << " MBps" << endl;
                }

                for (int i = 0; i < originalCount; ++i)
                {
                    blocks[i].Index = (unsigned char)i;
                    blocks[i].Block = originals[i];
                }

                for (int ii = 0; ii < recoveryCount; ++ii)
                {
                    int erasure_index = recoveryCount - ii - 1;
                    blocks[ii].Block = recoveryData + erasure_index * blockBytes;
                    blocks[ii].Index = (uint8_t)(originalCount + erasure_index);
                }

                {
                    LARGE_INTEGER t0; ::QueryPerformanceCounter(&t0);

                    if (cm256_decode(params, blocks))
                    {
                        cout << "Decoder error" << endl;
                        return 1;
                    }

                    LARGE_INTEGER t1; ::QueryPerformanceCounter(&t1);

                    LARGE_INTEGER tsum;
                    tsum.QuadPart = t1.QuadPart - t0.QuadPart;
                    double opusec = tsum.QuadPart * GetPerfFrequencyInverse() * 1000000.;
                    double mbps = (params.BlockBytes * params.OriginalCount / opusec);

                    cout << "Decoder: " << blockBytes << " bytes k = " << originalCount << " : " << opusec << " usec, " << mbps << " MBps" << endl;
                }

                if (!validateSolution(blocks, originalCount, blockBytes))
                {
                    cout << "Solution invalid" << endl;
                    return 1;
                }
            }
        }
    }

#if 0
    LARGE_INTEGER tsum;
    tsum.QuadPart = 0;

    cm256_encoder_params params;
    params.BlockBytes = 1296; // 1296
    params.OriginalCount = 100;
    params.RecoveryCount = 30;

    const int trials = 1000;
    for (int trial = 0; trial < trials; ++trial)
    {
        int x = 0;
        for (int i = 0; i < params.RecoveryCount * params.BlockBytes; ++i)
        {
            x ^= orig_data[i];
        }
        cout << x << endl;

        LARGE_INTEGER t0; ::QueryPerformanceCounter(&t0);

        if (cm256_encode(params, originals, recoveryData))
        {
            return 1;
        }

        LARGE_INTEGER t1; ::QueryPerformanceCounter(&t1);
        tsum.QuadPart += t1.QuadPart - t0.QuadPart;

        int y = 0;
        for (int i = 0; i < params.RecoveryCount * params.BlockBytes; ++i)
        {
            y ^= recoveryData[i];
        }
        cout << y << endl;
    }

    double opusec = tsum.QuadPart * GetPerfFrequencyInverse() * 1000000. / trials;
    double mbps = (params.BlockBytes * params.OriginalCount / opusec);

    cout << opusec << " usec, " << mbps << " MBps" << endl;
#endif

    return 0;
}
