## 1. Update usage text

- [x] 1.1 In `usage()`, rename `[-i HOST]` to `[-i AUDIO_HOST]` in the `--encode` usage line
- [x] 1.2 In `usage()`, rename `[-i HOST]` to `[-i AUDIO_HOST]` in the `--encode-image` usage line
- [x] 1.3 Add example line `./audio_steg --encode lsb "secret message" -i host.wav -o stego.wav` after the first encode example

## 2. Update file-header comment

- [x] 2.1 In the top-of-file comment block, rename `-i HOST.wav` to `-i AUDIO_HOST.wav` in the `--encode` mode description
- [x] 2.2 In the top-of-file comment block, rename `-i HOST.wav` to `-i AUDIO_HOST.wav` in the `--encode-image` mode description

## 3. Verify

- [x] 3.1 Build `audio_steg` and confirm the updated help output displays correctly
