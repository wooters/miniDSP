/**
 * @file test_fileio.c
 * @brief Tests for File I/O writers.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include "minidsp.h"
#include "fileio.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * Tests for FIO_write_npy()
 * -----------------------------------------------------------------------*/

/** Write a 3x4 float matrix to .npy and verify the raw bytes. */
static int test_write_npy(void)
{
    const char *fname = "_test_output.npy";
    unlink(fname);

    /* Build a known 3x4 matrix */
    float row0[] = {1.0f, 2.0f, 3.0f, 4.0f};
    float row1[] = {5.0f, 6.0f, 7.0f, 8.0f};
    float row2[] = {9.0f, 10.0f, 11.0f, 12.0f};
    const float *rows[] = {row0, row1, row2};

    if (FIO_write_npy(fname, (const float **)rows, 3, 4) != 0) {
        unlink(fname);
        return 0;
    }

    /* Read back the file */
    FILE *f = fopen(fname, "rb");
    if (!f) { unlink(fname); return 0; }

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    unsigned char *buf = malloc((size_t)fsize);
    fread(buf, 1, (size_t)fsize, f);
    fclose(f);

    int ok = 1;

    /* Check magic bytes */
    ok &= (buf[0] == 0x93);
    ok &= (buf[1] == 'N' && buf[2] == 'U' && buf[3] == 'M' &&
            buf[4] == 'P' && buf[5] == 'Y');
    /* Check version 1.0 */
    ok &= (buf[6] == 1 && buf[7] == 0);

    /* Read header length (LE u16) */
    uint16_t hlen = (uint16_t)(buf[8] | (buf[9] << 8));

    /* Check total prefix is divisible by 64 */
    ok &= ((10 + hlen) % 64 == 0);

    /* Check header contains shape string */
    char *header = (char *)(buf + 10);
    ok &= (strstr(header, "'shape': (3, 4)") != NULL);
    ok &= (strstr(header, "'descr': '<f4'") != NULL);

    /* Check raw float data */
    size_t data_offset = 10 + hlen;
    float *data = (float *)(buf + data_offset);
    float expected[] = {1,2,3,4, 5,6,7,8, 9,10,11,12};
    for (int i = 0; i < 12; i++) {
        ok &= approx_equal((double)data[i], (double)expected[i], 1e-6);
    }

    /* Check file size */
    ok &= ((size_t)fsize == data_offset + 12 * sizeof(float));

    free(buf);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for FIO_write_safetensors()
 * -----------------------------------------------------------------------*/

/** Write known data to safetensors and verify the raw bytes. */
static int test_write_safetensors(void)
{
    const char *fname = "_test_output.safetensors";
    unlink(fname);

    float row0[] = {1.0f, 2.0f, 3.0f};
    float row1[] = {4.0f, 5.0f, 6.0f};
    const float *rows[] = {row0, row1};

    if (FIO_write_safetensors(fname, (const float **)rows, 2, 3) != 0) {
        unlink(fname);
        return 0;
    }

    FILE *f = fopen(fname, "rb");
    if (!f) { unlink(fname); return 0; }

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    unsigned char *buf = malloc((size_t)fsize);
    fread(buf, 1, (size_t)fsize, f);
    fclose(f);

    int ok = 1;

    /* Read 8-byte LE u64 header size */
    uint64_t hsize = 0;
    for (int i = 0; i < 8; i++) {
        hsize |= ((uint64_t)buf[i]) << (i * 8);
    }

    /* JSON header should be reasonable size */
    ok &= (hsize > 10 && hsize < 500);

    /* Check JSON contains expected fields */
    char *json = malloc((size_t)hsize + 1);
    memcpy(json, buf + 8, (size_t)hsize);
    json[hsize] = '\0';

    ok &= (strstr(json, "\"dtype\":\"F32\"") != NULL);
    ok &= (strstr(json, "\"shape\":[2,3]") != NULL);
    ok &= (strstr(json, "\"data_offsets\":[0,24]") != NULL);  /* 2*3*4 = 24 */

    free(json);

    /* Check raw float data */
    size_t data_offset = 8 + (size_t)hsize;
    float *data = (float *)(buf + data_offset);
    float expected[] = {1,2,3, 4,5,6};
    for (int i = 0; i < 6; i++) {
        ok &= approx_equal((double)data[i], (double)expected[i], 1e-6);
    }

    /* Check file size */
    ok &= ((size_t)fsize == data_offset + 6 * sizeof(float));

    free(buf);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for FIO_write_wav()
 * -----------------------------------------------------------------------*/

/** Write a sine wave to WAV, read it back, verify round-trip. */
static int test_write_wav(void)
{
    const char *fname = "_test_output.wav";
    unlink(fname);

    unsigned samprate = 16000;
    size_t datalen = 1600;  /* 0.1 seconds */
    float *data = malloc(datalen * sizeof(float));

    /* Generate a 440 Hz sine wave */
    for (size_t i = 0; i < datalen; i++) {
        data[i] = (float)sin(2.0 * M_PI * 440.0 * (double)i / (double)samprate);
    }

    if (FIO_write_wav(fname, data, datalen, samprate) != 0) {
        free(data);
        unlink(fname);
        return 0;
    }

    /* Read it back */
    float *readback = NULL;
    size_t readlen = 0;
    unsigned readrate = 0;
    if (FIO_read_audio(fname, &readback, &readlen, &readrate, 0) != 0) {
        free(data);
        unlink(fname);
        return 0;
    }

    int ok = 1;

    /* Check sample count and rate */
    ok &= (readlen == datalen);
    ok &= (readrate == samprate);

    /* Check values round-trip within tolerance */
    if (readlen == datalen) {
        for (size_t i = 0; i < datalen; i++) {
            ok &= approx_equal((double)data[i], (double)readback[i], 1e-5);
        }
    }

    free(readback);
    free(data);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_fileio_tests(void)
{
    printf("\n--- File I/O writers ---\n");
    RUN_TEST(test_write_npy);
    RUN_TEST(test_write_safetensors);
    RUN_TEST(test_write_wav);
}
